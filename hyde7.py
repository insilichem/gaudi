###########################################################
# CHIMERA SCRIPT
# Rotate selected bonds, find clashes and H bonds
# By Jaime RG <jaime.rogue@gmail.com>, in UAB
###########################################################

## VERSION 7
# Implement genetic algorithm

# TODO
# - Consider Von Mises Distribution for random angle
# - Calculate h bonds only if clashes < threshold?
# - Avoid ligand rebuilding every step! Use 'ghost' library instead
# - Better crossover and mutations functions

# Chimera
import chimera, Rotamers, SwapRes
from chimera import UserError
# Python
import random, numpy, deap, argparse, sys, os
from deap import creator, tools, base, algorithms
# Custom
import hyde5, lego
import fragment3 as frag

### CUSTOM FUNCTIONS
def evalCoord(ind):

	## 1 - Build ligand
	lego.clearBase(cbase)
	linker = frag.insertMol(linkers[ind['molecule'][0]], target=anchor, join=True, 
		inplace=True, h=ind['h'][0])
	linker_anchor = [ a for a in linker if a.anchor in (4,6,8) ]
	frag.insertMol(fragments[int(ind['molecule'][1])], target=linker_anchor[0], 
		alpha=-120., h=ind['h'][1])		
	seen = { } # Used to remove duplicate entries in `bonds`
	bonds = [ seen.setdefault(b, b) for a in linker for b in a.bonds if b not in seen ]

	## 2 - Set rotations
	# Direct access to BondRot, instead of BondRotMgr
	for i, bond in enumerate(bonds):
		hyde5.bondrot(bond, anchor=cbase[0], delta=ind['linker_rots'][i])

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
	all_atoms = [ a for m in chimera.openModels.list() for a in m.atoms ]
	res_atoms = [ a for r in residues for a in r.atoms ]
	ligand_atoms = cbase[0].residue.atoms
	not_ligand = list(set(all_atoms) - set(ligand_atoms))

	# TODO: Restrict donor and acceptors to smaller selection
	hbonds = hyde5.countHBonds(model, sel=ligand_atoms, cache=False)
	clashes, num_of_clashes = hyde5.countClashes(atoms=ligand_atoms, 
		test=all_atoms) # TODO: Restrict test to smaller selection
	clashes_r, num_of_clashes_r = hyde5.countClashes(atoms=res_atoms,
		test=not_ligand)

	return len(hbonds), num_of_clashes, num_of_clashes_r

def hetCxOnePoint(ind1, ind2):

	for key in ind1:
		if key == 'molecule': #ignore building blocks
			continue
		size = min(len(ind1[key]), len(ind2[key]))
		if size > 1:
			cxpoint = random.randint(1, size - 1)
			ind1[key][cxpoint:], ind2[key][cxpoint:] = \
			ind2[key][cxpoint:], ind1[key][cxpoint:]

	return ind1, ind2

def hetMutation(ind, indpb):
	for key, row in ind.items():
		if random.random() < ind:
			j = random.randint(0,len(row)-1)
			if key == 'molecule': 
				continue
			elif key == 'h':
				ind[key][j] = toolbox.rand_h()	
			elif key == 'linker_rots':
				ind[key][j] = toolbox.rand_angle()
			elif key == 'mutamers':
				ind[key][j] = toolbox.rand_aa()
			elif key == 'rotamers':
				ind[key][j] = toolbox.rand_rotamer()		
	return ind,

aminoacids = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLU', 'GLN',
			  'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE',
			  'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL' ]
##/ FUNCTIONS

# Initial checks
required_sels = 'base', 'anchor', 'rotamers'
if not set(required_sels) <= set(chimera.selection.savedSels):
	raise UserError("""You have to define three selections with 'namesel'
'base': Single atom. Terminal end of the static part of the ligand
'anchor': Single atom. Together with 'base', it will define the static part.
	It's also the insertion point for the linker.
'mutable': Residues to be swapped and mutated.
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
args = parser.parse_args()
#### /ARGUMENT PARSING

wd = os.path.dirname(os.path.realpath(sys.argv[0]))
# Get Chimera params
model = chimera.openModels.list()
residues = chimera.selection.savedSels['mutable'].residues()
base_at = chimera.selection.savedSels['base'].atoms()[0]
anchor = chimera.selection.savedSels['anchor'].atoms()[0]

# Get building blocks
cbase = [base_at] + list(hyde5.atomsBetween(base_at, anchor)) + [anchor]
linkers = lego.getMol2Files(wd +'/mol2/linkers/')
fragments = lego.getMol2Files(wd + '/mol2/fragments/')

###
# Genetic Algorithm
# define individual, population, etc
deap.creator.create("FitnessMax", deap.base.Fitness, weights=(1.0, -1.0, -1.0))
deap.creator.create("Individual", dict, fitness=deap.creator.FitnessMax)

# Operators
toolbox = deap.base.Toolbox()
toolbox.register("rand_angle", random.uniform, 0, 360)
toolbox.register("rand_h", random.randint, 1, 3)
toolbox.register("rand_aa", random.randint, 0, len(aminoacids)-1)
toolbox.register("rand_rotamer", random.randint, 0, 8)
toolbox.register("rand_linker", random.randint, 0, len(linkers)-1)
toolbox.register("rand_fragment", random.randint, 0, len(fragments)-1)

# Genes
toolbox.register("molecule", deap.tools.initCycle, list,
	[toolbox.rand_linker, toolbox.rand_fragment], n=1)
toolbox.register("h", deap.tools.initRepeat, list,
	toolbox.rand_h, n=2)
toolbox.register("linker_rots", deap.tools.initRepeat, list,
	toolbox.rand_angle, n=8)
# toolbox.register("fragment_rots", deap.tools.initRepeat, list,
# 	toolbox.rand_angle)
toolbox.register("mutamers", deap.tools.initRepeat, list,
	toolbox.rand_aa, n=len(residues))
toolbox.register("rotamers", deap.tools.initRepeat, list,
	toolbox.rand_rotamer, n=len(residues))
toolbox.register("toDict", 
	(lambda ind, *fn: ind((f.__name__, f()) for f in fn)))

# Individual and population
toolbox.register("individual", toolbox.toDict, deap.creator.Individual, 
	toolbox.molecule, toolbox.h, toolbox.linker_rots,
	toolbox.mutamers, toolbox.rotamers)
toolbox.register("population", deap.tools.initRepeat, 
	list, toolbox.individual)

# Aliases for algorithm
toolbox.register("evaluate", evalCoord)
toolbox.register("mate", hetCxOnePoint)
toolbox.register("mutate", hetMutation, indpb=0.05)
toolbox.register("select", deap.tools.selNSGA2)

def main():
	pop = toolbox.population(n=args.pop-1)
	hof = deap.tools.HallOfFame(1)
	stats = deap.tools.Statistics(lambda ind: ind.fitness.values[1])
	stats.register("avg", numpy.mean)
	stats.register("std", numpy.std)
	stats.register("min", numpy.min)
	stats.register("max", numpy.max)
	pop, log = deap.algorithms.eaMuPlusLambda(pop, toolbox, 
		mu = int(args.pop/2), lambda_= int(args.pop/2), cxpb=0.5, 
		mutpb=0.2, ngen=args.ngen, stats=stats, halloffame=hof)
	return pop, log, hof

if __name__ == "__main__":
	pop, log, hof = main()
	evalCoord(hof[0])
	print("Best individual is: %s\nwith fitness: %s" % (hof[0], hof[0].fitness))
	# test = {'h': [3, 3], 
	# 	'rotamers': [4, 5], 
	# 	'mutamers': [12, 9], 
	# 	'linker_rots': [268, 44, 207, 298, 248, 355, 200, 297], 
	# 	'molecule': [2, 0]}
	# print "Individual:\n{0}\nFitness:\n{1}".format(test, evalCoord(test))
