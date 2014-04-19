###########################################################
# CHIMERA SCRIPT
# Rotate selected bonds, find clashes and H bonds
# By Jaime RG <jaime.rogue@gmail.com>, in UAB
###########################################################

## VERSION 7
# Implement genetic algorithm

# Chimera
import chimera, Rotamers, SwapRes
from chimera import UserError
# Python
import random, numpy, deap, sys
from deap import creator, tools, base, algorithms
# Custom
import mof3d
from mof3d.utils import box
reload(mof3d)
### CUSTOM FUNCTIONS

def evalCoord(ind, close=True, hidden=False):

	## 1 - Choose ligand from pre-built mol library
	ligand, bondrots = ligands[ind['ligand']]
	chimera.openModels.add([ligand], shareXform=True, hidden=hidden)
	# Chimera converts metal bonds to pseudoBonds all the time
	box.pseudobond_to_bond(ligand)

	if 'rotable_bonds' in ind:
		for alpha, br in zip(ind['rotable_bonds'], bondrots):
			chimera.selection.addCurrent([br.bond, br.rotanchor])
			br.adjustAngle(alpha - br.angle, br.rotanchor)
	
	if 'rotamers' in ind:
		if 'mutamers' in ind:
			aas = [ AA[i] for i in ind['mutamers'] ]
		else:
			aas = [ r.type for r in residues ]
		
		for res, rot, mut in map(None, residues, ind['rotamers'], aas):
			try: 
				rotamers = Rotamers.getRotamers(res, resType=mut, 
												lib=cfg.rotamers.library.title())[1]
				Rotamers.useRotamer(residues[i],[rotamers[rot]])
			except Rotamers.NoResidueRotamersError: # ALA, GLY...
				if 'mutamers' in ind:
					SwapRes.swap(residues[i], AA[mut], bfactor=None)
			except IndexError:
				Rotamers.useRotamer(residues[i],rotamers[-1:])

	ligand_env.clear()
	ligand_env.add(ligand)
	ligand_env.merge(chimera.selection.REPLACE, 
					chimera.specifier.zone( ligand_env, 'atom', None, 15.0, 
											[protein,ligand]))
	score = []
	for obj in cfg.objective:
		if obj.type == 'hbonds':
			hbonds = mof3d.score.chem.hbonds(
						[protein, ligand], cache=False, test=ligand_env.atoms(),
						sel=[ a for a in ligand.atoms if a not in ("C", "CA", "N", "O") ])

			score.append(len(hbonds))
		
		elif obj.type in ('clashes', 'contacts'):
			contacts, num_of_contacts, positive_vdw, negative_vdw =\
				mof3d.score.chem.clashes(atoms=ligand.atoms, 
										test=ligand_env.atoms(), 
										intraRes=True, clashThreshold=obj.threshold, 
										hbondAllowance=0.0, parse=True)
			if obj.type == 'clashes':
				score.append(sum(abs(a[3]) for a in negative_vdw)/2)
			else:
				score.append(sum(1-a[3] for a in positive_vdw)/2)
			
			if not close:
				if obj.type == 'contacts' and positive_vdw: 
					mof3d.score.chem.draw_clashes(positive_vdw, startCol=obj.color[0],
						endCol=obj.color[1], key=3, name="Hydrophobic interactions")
				if obj.type == 'clashes' and negative_vdw:
					mof3d.score.chem.draw_clashes(negative_vdw, startCol=obj.color[0], 
						endCol=obj.color[1], key=3, name="Clashes")
		
		elif obj.type == 'distance' :
			probes = box.atoms_by_serial(*obj.probes, atoms=ligand.atoms)
			target, = box.atoms_by_serial(obj.target, atoms=protein.atoms)
			dist = mof3d.score.target.distance(probes, target, obj.threshold, \
				wall=obj.wall)
			score.append(dist)
	if close:
		chimera.openModels.remove([ligand])

	return score

def hetCrossover(ind1, ind2):
	for key in ind1:
		if key == 'molecule': 
			continue
		elif key == 'rotable_bonds':
			ind1[key][:], ind2[key][:] = deap.tools.cxSimulatedBinaryBounded(
				ind1[key], ind2[key], eta=10., low=0., up=360.)
		elif key == 'mutamers':
			ind1[key], ind2[key] = deap.tools.cxTwoPoint(ind1[key], ind2[key])
		elif key == 'rotamers':
			ind1[key], ind2[key] = deap.tools.cxTwoPoint(ind1[key], ind2[key])
	
	return ind1, ind2

def hetMutation(ind, indpb,):
	for key in ind:
		if key == 'molecule': 
			continue
		elif key == 'rotable_bonds':
			ind[key] = deap.tools.mutPolynomialBounded(ind[key], 
				eta=10., low=0., up=360., indpb=indpb)[0]
		elif key == 'mutamers':
			ind[key] = deap.tools.mutUniformInt(ind[key], 
				low=0, up=len(residues)-1, indpb=indpb)[0]
		elif key == 'rotamers':
			ind[key] = deap.tools.mutUniformInt(ind[key], 
				low=0, up=8, indpb=indpb)[0]
	
	return ind,

##/ FUNCTIONS

## Initialize workspace
cfg = mof3d.utils.parse.Settings(sys.argv[1])
deap.creator.create("FitnessMax", deap.base.Fitness, weights=cfg.weights())
deap.creator.create("Individual", dict, fitness=deap.creator.FitnessMax)

protein, = chimera.openModels.open(cfg.protein.path)
ligand_env = chimera.selection.ItemizedSelection()

if cfg.ligand.type == 'mol2':
	ligand = chimera.openModels.open(cfg.ligand.path)
	ligands = { id(ligand): ligand }
elif cfg.ligand.type == 'blocks':
	rotations = True if cfg.ligand.flexible else False
	ligands = mof3d.molecule.library(cfg.ligand.path, 
				bondto=box.atoms_by_serial(cfg.ligand.bondto, atoms=protein.atoms)[0],
				rotations=rotations, join='dummy')

# Operators and genes
genes = []
toolbox = deap.base.Toolbox()
toolbox.register("ligand", random.choice, ligands.keys())
genes.append(toolbox.ligand)

if hasattr(cfg.ligand, 'rotable') or (hasattr(cfg.ligand, 'bondto') and cfg.ligand.bondto):
	toolbox.register("rand_angle", random.uniform, 0, 360)
	toolbox.register("rotable_bonds", deap.tools.initRepeat, list,
						toolbox.rand_angle, n=20)
	genes.append(toolbox.rotable_bonds) 

if not hasattr(cfg.ligand, 'bondto') or not cfg.ligand.bondto:
	pass #activate 3D global search

if hasattr(cfg, 'rotamers'):
	residues = [ r for r in protein.residues if r.id.position in cfg.rotamers.residues ]
	toolbox.register("rand_rotamer", random.randint, 0, cfg.rotamers.top-1)
	toolbox.register("rotamers", deap.tools.initRepeat, list,
					toolbox.rand_rotamer, n=len(cfg.rotamers.residues))
	genes.append(toolbox.rotamers)
	if cfg.rotamers.mutate == "all" or isinstance(cfg.rotamers.mutate, list):
		AA = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLU', 'GLN',
			  'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE',
			  'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL' ]
		if isinstance(cfg.rotamers.mutate, list):
			AA = [ AA[i-1] for i in cfg.rotamers.mutate if i<=20 ]
		toolbox.register("rand_aa", random.randint, 0, len(AA)-1)
		toolbox.register("mutamers", deap.tools.initRepeat, list,
						toolbox.rand_aa, n=len(cfg.rotamers.residues))
		genes.append(toolbox.mutamers)

toolbox.register("toDict", 
	(lambda ind, *fn: ind((f.__name__, f()) for f in fn)))

# Individual and population
toolbox.register("individual", toolbox.toDict, deap.creator.Individual, *genes)
toolbox.register("population", deap.tools.initRepeat, list, toolbox.individual)

# Aliases for algorithm
toolbox.register("evaluate", evalCoord)
toolbox.register("mate", hetCrossover)
toolbox.register("mutate", hetMutation, indpb=0.05)
toolbox.register("select", deap.tools.selNSGA2)


def main():
	pop = toolbox.population(n=cfg.ga.pop)
	hof = deap.tools.ParetoFront() if cfg.ga.pareto \
			else deap.tools.HallOfFame(cfg.default.results)
	stats = deap.tools.Statistics(lambda ind: ind.fitness.values)
	numpy.set_printoptions(precision=cfg.default.precision)
	stats.register("avg", numpy.mean, axis=0)
	stats.register("min", numpy.min, axis=0)
	stats.register("max", numpy.max, axis=0)
	pop, log = deap.algorithms.eaMuPlusLambda(pop, toolbox, 
		mu = cfg.ga.mu, lambda_= cfg.ga.lambda_, cxpb=cfg.ga.cx_pb, 
		mutpb=cfg.ga.mut_pb, ngen=cfg.ga.gens, stats=stats, halloffame=hof)
	return pop, log, hof

if __name__ == "__main__":	
	print "Scores: " + ', '.join([o.type for o in cfg.objective])
	pop, log, hof = main()
	evalCoord(hof[0], close=False)
	print("Best individual is: %s\nwith fitness: %s" % (hof[0], hof[0].fitness))
	print("More possible solutions to assess:")
	for h in hof[1:]:
		print h, h.fitness
	# for i, (ligand, br) in ligands.items():
	# 	chimera.openModels.add([ligand])