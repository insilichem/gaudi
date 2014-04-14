###########################################################
# CHIMERA SCRIPT
# Rotate selected bonds, find clashes and H bonds
# By Jaime RG <jaime.rogue@gmail.com>, in UAB
###########################################################

## VERSION 7
# Implement genetic algorithm

# TODO
# - Calculate h bonds only if clashes < threshold?


# Chimera
import chimera, Rotamers, SwapRes, SplitMolecule, ChemGroup as cg
from chimera import UserError
# Python
import random, numpy, deap, argparse, sys, os
from deap import creator, tools, base, algorithms
# Custom
import fetra
### CUSTOM FUNCTIONS

def evalCoord(ind, close=True, hidden=False):
	## 1 - Choose ligand from pre-built mol library
	ligand, bondrots = mol_library[ind['molecule'][0],ind['molecule'][1]]
	chimera.openModels.add([ligand], shareXform=True, hidden=hidden)
	# Chimera converts metal bonds to pseudoBonds all the time
	fetra.utils.box.pseudobond_to_bond(ligand)

	## 2 - Set rotations
	# Direct access to BondRot, instead of BondRotMgr
	for i, br in enumerate(bondrots):
		br.adjustAngle(ind['linker_rots'][i] - br.angle, br.myanchor)

	## 3 - Set mutamers/rotamers
	for i, aa in enumerate(ind['mutamers']):
		try: 
			rotamers = Rotamers.getRotamers(residues[i], resType=aminoacids[aa])[1] 
			rotId = ind['rotamers'][i]
			Rotamers.useRotamer(residues[i],[rotamers[rotId]])
		except Rotamers.NoResidueRotamersError: # ALA, GLY...
			SwapRes.swap(residues[i], aminoacids[aa], bfactor=None)
		except IndexError:
			Rotamers.useRotamer(residues[i],[rotamers[-1]])

	## 4 - Score
	# Inserted atoms become new entities
	res_atoms = [ a for r in residues for a in r.atoms ]
	models = chimera.openModels.list(all=True,modelTypes=[chimera.Molecule])
	ligand_env.clear()
	ligand_env.add(ligand)
	ligand_env.merge(chimera.selection.REPLACE, 
					chimera.specifier.zone(ligand_env, 'atom', None, 15.0, models))

	hbonds = fetra.score.chem.hbonds(
				models, cache=False, test=ligand_env.atoms(),
				sel=[ a for a in ligand.atoms if a not in ("C", "CA", "N", "O") ])
	# TODO: Restrict test to smaller selection
	contacts, num_of_contacts, positive_vdw, negative_vdw =\
		fetra.score.chem.clashes(atoms=ligand.atoms, 
								test=ligand_env.atoms(), 
								intraRes=True, clashThreshold=-0.4, 
								hbondAllowance=0.0, parse=True)
	clashes_r, num_of_clashes_r = fetra.score.chem.clashes(atoms=res_atoms,
														test=mol.atoms)
	
	if not close:
		if positive_vdw: 
			fetra.score.chem.draw_clashes(positive_vdw, startCol="FFFF00", endCol="00FF00",
				key=3, name="Hydrophobic interactions")
		if negative_vdw:
			fetra.score.chem.draw_clashes(negative_vdw, startCol="FF0000", endCol="FF0000",
				key=3, name="Bad clashes")
	else:
		chimera.openModels.remove([ligand])

	return  len(hbonds), sum(abs(a[3]) for a in negative_vdw)/2, \
			sum(1-a[3] for a in positive_vdw)/2, num_of_clashes_r

def hetCrossover(ind1, ind2):
	for key in ind1:
		if key == 'molecule': 
			continue
		elif key == 'linker_rots':
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
		elif key == 'linker_rots':
			ind[key] = deap.tools.mutPolynomialBounded(ind[key], 
				eta=10., low=0., up=360., indpb=indpb)[0]
		elif key == 'mutamers':
			ind[key] = deap.tools.mutUniformInt(ind[key], 
				low=0, up=len(residues)-1, indpb=indpb)[0]
		elif key == 'rotamers':
			ind[key] = deap.tools.mutUniformInt(ind[key], 
				low=0, up=8, indpb=indpb)[0]
	return ind,

# Couple of constants
aminoacids = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLU', 'GLN',
			  'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE',
			  'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL' ]
##/ FUNCTIONS

# Initial checks
required_sels = 'base', 'anchor', 'mutable'
if not set(required_sels) <= set(chimera.selection.savedSels):
	raise UserError("""You have to define three selections with 'namesel'
---
- 'base': Single atom. Terminal end of the static part of the ligand
- 'anchor': Single atom. Together with 'base', it will define the static part. It's also the insertion point for the linker.
- 'mutable': Residues to be swapped and mutated.
""")

parser = argparse.ArgumentParser()
parser.add_argument('-p', '--population',
					required=False,
					type=int,
					default=100, 
					dest="pop",
					metavar="<number of individuals>",
					help="Number of individuals in each population" )
parser.add_argument('-g', '--generation',
					required=False,
					type=int,
					default=30, 
					dest="ngen",
					metavar="<number of generations>",
					help="Number of generations to calculate" )
parser.add_argument('-e', '--eval',
					required=False,
					type=str,
					default='', 
					dest="eval",
					metavar="<individual>",
					help="Paste individual data to evaluate" )
parser.add_argument('-pf', '--pareto',
					required=False,
					default=False, 
					dest="pareto",
					action="store_true",
					help="Wether to report Pareto Front rank or not" )
args = parser.parse_args()
wd = os.path.dirname(os.path.realpath(sys.argv[0]))

# Get Chimera params
residues = [ r for r in sorted(chimera.selection.savedSels['mutable'].residues(), 
	key= lambda r: (r.id.chainId, r.id.position)) ]
base_at = chimera.selection.savedSels['base'].atoms()[0]
mol = base_at.molecule
ligand = base_at.residue
anchor = chimera.selection.savedSels['anchor'].atoms()[0]
cbase = [base_at] + list(fetra.utils.box.atoms_between(base_at, anchor)) + [anchor]
ligand_env = chimera.selection.ItemizedSelection()
#dihedral
d3 = anchor
d4 = [a for a in anchor.neighbors if a not in cbase][0]
d2 = [a for a in d3.neighbors if a in cbase][0]
d1 = [a for a in d2.neighbors if a != anchor ]
if len(d1)>1: 
	d1 = [a for a in d1 if a.numBonds>1]
d1 = d1[0]
dihedral = chimera.dihedral(d1.coord(),d2.coord(),d3.coord(),d4.coord())
#angle
alpha = chimera.angle(d2.coord(),d3.coord(),d4.coord())

# Save cbase in ghost molecule
cbasecopy = SplitMolecule.split.molecule_from_atoms(mol, cbase)
[mol.deleteAtom(a) for a in ligand.atoms]

# Get building blocks
linkers = sorted(fetra.utils.box.files_in(wd +'/mol2/linkers/', 'mol2'))
fragments = sorted(fetra.utils.box.files_in(wd +'/mol2/fragments/', 'mol2'))

# Build library
mol_library = fetra.molecule.library(cbasecopy, linkers, fragments, \
									dihedral=dihedral, alpha=alpha)

###
# Genetic Algorithm
# define individual, population, etc
deap.creator.create("FitnessMax", deap.base.Fitness, weights=(1.0, -1.0, 1.0, -1.0))
deap.creator.create("Individual", dict, fitness=deap.creator.FitnessMax)

# Operators
toolbox = deap.base.Toolbox()
toolbox.register("rand_angle", random.uniform, 0, 360)
toolbox.register("rand_aa", random.randint, 0, len(aminoacids)-1)
toolbox.register("rand_rotamer", random.randint, 0, 8)
toolbox.register("rand_linker", random.randint, 0, len(linkers)-1)
toolbox.register("rand_fragment", random.randint, 0, len(fragments)-1)

# Genes
toolbox.register("molecule", deap.tools.initCycle, list,
	[toolbox.rand_linker, toolbox.rand_fragment], n=1)
toolbox.register("linker_rots", deap.tools.initRepeat, list,
	toolbox.rand_angle, n=8)
toolbox.register("mutamers", deap.tools.initRepeat, list,
	toolbox.rand_aa, n=len(residues))
toolbox.register("rotamers", deap.tools.initRepeat, list,
	toolbox.rand_rotamer, n=len(residues))
toolbox.register("toDict", 
	(lambda ind, *fn: ind((f.__name__, f()) for f in fn)))

# Individual and population
toolbox.register("individual", toolbox.toDict, deap.creator.Individual, 
	toolbox.molecule, toolbox.linker_rots, toolbox.mutamers, toolbox.rotamers)
toolbox.register("population", deap.tools.initRepeat, list, toolbox.individual)

# Aliases for algorithm
toolbox.register("evaluate", evalCoord)
toolbox.register("mate", hetCrossover)
toolbox.register("mutate", hetMutation, indpb=0.05)
toolbox.register("select", deap.tools.selNSGA2)

def main():
	pop = toolbox.population(n=args.pop)
	if args.pareto:
		hof = deap.tools.ParetoFront()
	else: 
		hof = deap.tools.HallOfFame(20)
	stats = deap.tools.Statistics(lambda ind: ind.fitness.values)
	stats.register("avg", numpy.mean, axis=0)
	stats.register("min", numpy.min, axis=0)
	stats.register("max", numpy.max, axis=0)
	pop, log = deap.algorithms.eaMuPlusLambda(pop, toolbox, 
		mu = int(args.pop/2), lambda_= int(args.pop/2), cxpb=0.5, 
		mutpb=0.2, ngen=args.ngen, stats=stats, halloffame=hof)
	return pop, log, hof

if __name__ == "__main__":
	numpy.set_printoptions(precision=2)
	if args.eval:
		print "Fitness: ", evalCoord(eval(args.eval), close=False)
	else:
		pop, log, hof = main()
		evalCoord(hof[0], close=False)
		print("Best individual is: %s\nwith fitness: %s" % (hof[0], hof[0].fitness))
		print("More possible solutions to assess:")
		for h in hof[1:11]:
			print h, h.fitness
