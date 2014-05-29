###########################################################
# CHIMERA SCRIPT
# Rotate selected bonds, find clashes and H bonds
# By Jaime RG <jaime.rogue@gmail.com>, in UAB
###########################################################

## VERSION 7
# Implement genetic algorithm

# Chimera
import chimera, Rotamers, SwapRes, Matrix as M
# Python
import random, numpy, deap, sys, math
from deap import creator, tools, base, algorithms
from operator import mul
# Custom
import gaudi
from gaudi.utils import box

### CUSTOM FUNCTIONS

def evaluate(ind, close=True, hidden=False, draw=False):

	ligand, bondrots = ligands[ind['ligand']]
	chimera.openModels.add([ligand], shareXform=True, hidden=hidden)
	box.pseudobond_to_bond(ligand)
	if 'rotable_bonds' in ind:
		for alpha, br in zip(ind['rotable_bonds'], bondrots):
			# check amides!
			if all([a.idatmType in ('C2', 'N2') for a in br.bond.atoms]):
				alpha = 0. if alpha<180 else 180.
			br.adjustAngle(alpha - br.angle, br.rotanchor)

	if 'xform' in ind:
		ligand.openState.xform = M.chimera_xform(M.multiply_matrices(*ind['xform']))

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
					SwapRes.swap(res, mut, bfactor=None)
			except IndexError: #use last available
				Rotamers.useRotamer(res,rotamers[-1:])

	ligand_env.clear()
	ligand_env.add(ligand.atoms)
	ligand_env.merge(chimera.selection.REPLACE,
					chimera.specifier.zone( ligand_env, 'atom', None, 10.0,
											[protein,ligand]))
	score, draw_list = [], {}
	for obj in cfg.objective:
		if obj.type == 'hbonds':
			hbonds = gaudi.score.chem.hbonds(
						[protein, ligand], cache=False, test=ligand_env.atoms(),
						sel=[ a for a in ligand.atoms if a not in ("C", "CA", "N", "O") ],
						dist_slop=obj.distance_tolerance, angle_slop=obj.angle_tolerance)
			score.append(len(hbonds))
			draw_list['hbonds'] = hbonds
			if hasattr(obj, 'targets') and len(obj.targets):
				hb_targets = box.atoms_by_serial(*obj.targets, atoms=protein.atoms)
				hb_targets_num = [ha for hb in hbonds for ha in hb if ha in hb_targets]
				score.append(len(hb_targets_num))
		
		elif obj.type == 'contacts':
			if cfg.ligand.type == 'blocks':
				ligand_atoms = [a for a in ligand.atoms if a.serialNumber > 3]
			else:
				ligand_atoms = ligand.atoms
			contacts, num_of_contacts, positive_vdw, negative_vdw = \
				gaudi.score.chem.clashes(atoms=ligand_atoms, 
										test=ligand_env.atoms(), 
										intraRes=True, clashThreshold=obj.threshold, 
										hbondAllowance=obj.threshold_h, parse=True,
										parse_threshold=obj.threshold_c)
			if obj.which == 'clashes':
				clashscore = sum(abs(a[3]) for a in negative_vdw)/2
				if clashscore > obj.cutoff and close:
					chimera.openModels.remove([ligand])
					return [-1000*w for w in weights]
				score.append(clashscore)
				draw_list['negvdw'] = negative_vdw
			elif obj.which == 'hydrophobic':
				score.append(sum(1-a[3] for a in positive_vdw)/2)
				draw_list['posvdw'] = positive_vdw

		elif obj.type == 'distance':
			probes = box.atoms_by_serial(*obj.probes, atoms=ligand.atoms)
			dist_target, = box.atoms_by_serial(obj.target, atoms=protein.atoms)
			dist = gaudi.score.target.distance(probes, dist_target, obj.threshold, 
						wall=obj.wall)
			score.append(dist)

		elif obj.type == 'solvation':
			atoms, ses, sas = gaudi.score.chem.solvation(ligand_env.atoms())
			ligand_ses = sum(s for (a,s) in zip(atoms, ses) if a in ligand.atoms)
			score.append(ligand_ses)

	if close:
		chimera.openModels.remove([ligand])
		return score
	if draw and draw_list:
		if 'negvdw' in draw_list:
			gaudi.score.chem.draw_interactions(draw_list['negvdw'], startCol='FF0000', 
				endCol='FF0000', key=3, name="Clashes")
		if 'posvdw' in draw_list: 
			gaudi.score.chem.draw_interactions(draw_list['posvdw'], startCol='00FF00',
				endCol='FFFF00', key=3, name="Hydrophobic interactions")
		if 'hbonds' in draw_list:
			gaudi.score.chem.draw_interactions(draw_list['hbonds'], startCol='00FFFF', 
				endCol='00FFFF', name="H Bonds")
	return ligand

def het_crossover(ind1, ind2):
	for key in ind1:
		if key == 'rotable_bonds':
			ind1[key][:], ind2[key][:] = deap.tools.cxSimulatedBinaryBounded(
				ind1[key], ind2[key], eta=cfg.ga.cx_eta, 
				low=-0.5*cfg.ligand.flexibility, up=0.5*cfg.ligand.flexibility)
		elif key in ('mutamers', 'rotamers'):
			ind1[key], ind2[key] = deap.tools.cxTwoPoint(ind1[key], ind2[key])
		elif key == 'xform': # swap rotation and translation
			ind1[key], ind2[key] = ind1[key][:1]+ind2[key][1:], ind2[key][:1]+ind1[key][1:]
			# ind1[key], ind2[key] = gaudi.move.pos_swap(ind1[key], ind2[key])

	return ind1, ind2

def het_mutation(ind, indpb):
	for key in ind:
		if key == 'rotable_bonds':
			ind[key] = deap.tools.mutPolynomialBounded(ind[key], 
				eta=cfg.ga.mut_eta, low=-0.5*cfg.ligand.flexibility, 
				up=0.5*cfg.ligand.flexibility, indpb=indpb)[0]
		elif key == 'mutamers':
			ind[key] = deap.tools.mutUniformInt(ind[key], 
				low=0, up=len(residues)-1, indpb=indpb)[0]
		elif key == 'rotamers':
			ind[key] = deap.tools.mutUniformInt(ind[key], 
				low=0, up=8, indpb=indpb)[0]
		elif key == 'xform' and random.random() < indpb:
			# Careful! Mutation generates a whole NEW position (similar to eta ~= 0)
			ind['xform'] = gaudi.move.rand_xform(origin, cfg.protein.radius)

	return ind,

def similarity(a, b):
	atoms1, atoms2 = ligands[a['ligand']][0].atoms, ligands[b['ligand']][0].atoms
	atoms1.sort(key=lambda x: x.serialNumber)
	atoms2.sort(key=lambda x: x.serialNumber)

	xf1, xf2 = M.multiply_matrices(*a['xform']), M.multiply_matrices(*b['xform'])
	xf1, xf2 = M.chimera_xform(xf1), M.chimera_xform(xf2)
	
	sqdist = sum( xf1.apply(a.coord()).sqdistance(xf2.apply(a.coord())) 
					for a, b in zip(atoms1, atoms2) )
	rmsd = math.sqrt(sqdist / float(len(atoms1)))
	
	return rmsd < 0.25
##/ FUNCTIONS

## Initialize workspace
cfg = gaudi.utils.parse.Settings(sys.argv[1])
weights = cfg.weights() if len(sys.argv)<=2 else map(float, sys.argv[2:])
deap.creator.create("FitnessMax", deap.base.Fitness, weights=(1.0,))
deap.creator.create("Individual", dict, fitness=deap.creator.FitnessMax,
					fitness_names=[o.type for o in cfg.objective])

# Open protein
protein, = chimera.openModels.open(cfg.protein.path)

# Set up ligands
ligand_env = chimera.selection.ItemizedSelection()
if not hasattr(cfg.ligand, 'bondto') or not cfg.ligand.bondto:
	target = origin = box.atoms_by_serial(cfg.protein.origin, atoms=protein.atoms)[0].coord()
	search3D = True
	join = False

else:
	target = box.atoms_by_serial(cfg.ligand.bondto, atoms=protein.atoms)[0]
	search3D = False
	join = 'dummy'
	
rotations = True if cfg.ligand.flexibility else False
ligands = gaudi.molecule.library(cfg.ligand.path,
				bondto=target, rotations=rotations, join=join)

# Operators and genes
genes = []
toolbox = deap.base.Toolbox()
toolbox.register("ligand", random.choice, ligands.keys())
genes.append(toolbox.ligand)

if search3D:
	toolbox.register("xform", gaudi.move.rand_xform, origin, cfg.protein.radius)
	genes.append(toolbox.xform)

if hasattr(cfg.ligand, 'flexibility') and cfg.ligand.flexibility:
	if cfg.ligand.flexibility > 360: 
		cfg.ligand.flexibility = 360.0
	toolbox.register("rand_angle", random.uniform, -0.5*cfg.ligand.flexibility, 0.5*cfg.ligand.flexibility)
	toolbox.register("rotable_bonds", deap.tools.initRepeat, list,
						toolbox.rand_angle, n=20)
	genes.append(toolbox.rotable_bonds) 

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
toolbox.register("evaluate", evaluate)
toolbox.register("mate", het_crossover)
toolbox.register("mutate", het_mutation, indpb=cfg.ga.mut_indpb)
toolbox.register("select", deap.tools.selNSGA2)

def main():
	pop = toolbox.population(n=cfg.ga.pop)
	hof = deap.tools.ParetoFront(similarity) if cfg.ga.pareto \
			else deap.tools.HallOfFame(cfg.default.results, similarity)
	stats = deap.tools.Statistics(lambda ind: ind.fitness.values)
	numpy.set_printoptions(precision=cfg.default.precision)
	stats.register("avg", numpy.mean, axis=0)
	stats.register("min", numpy.min, axis=0)
	stats.register("max", numpy.max, axis=0)
	pop, log = deap.algorithms.eaMuPlusLambda(pop, toolbox, 
		mu = int(cfg.ga.mu*cfg.ga.pop), lambda_= int(cfg.ga.lambda_*cfg.ga.pop), 
		cxpb=cfg.ga.cx_pb, mutpb=cfg.ga.mut_pb, 
		ngen=cfg.ga.gens, stats=stats, halloffame=hof)
	return pop, log, hof

if __name__ == "__main__":	
	print "Scores:", ', '.join(o.type for o in cfg.objective)
	pop, log, hof = main()

	rank = box.write_individuals(hof, cfg.default.savepath,	cfg.default.savename, evaluate)
	with open(cfg.default.savepath+cfg.default.savename+'.gaudi', 'w+') as out:
		out.write('\n'.join(rank))

	#Display best individual
	if not chimera.nogui:
		evaluate(hof[0], close=False, draw=True)