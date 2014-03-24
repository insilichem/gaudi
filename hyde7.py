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

import chimera, Rotamers, random, numpy, deap, argparse, sys, os
from deap import creator, tools, base, algorithms
from chimera import UserError
import hyde5, lego # the scripts!
import fragment3 as frag
reload(frag)
reload(hyde5)
### CUSTOM FUNCTIONS
#
def evalCoord(ind):

	## 1 - Build ligand
	# 'base' has to be selected in Chimera
	# 'anchor' should be generated on the fly
	cbase, anchor = lego.getBase()
	lego.clearBase(cbase, anchor)

	linker = frag.insertMol(linkers[int(ind[0][0])], target=anchor, join=True, 
		inplace=True, h=int(ind[1][0]))
	linker_anchor = [ a for a in linker if a.anchor in (4,6,8) ]
	fragm = frag.insertMol(fragments[int(ind[0][1])], target=linker_anchor[0], 
		alpha=-120., h=int(ind[1][1]))		
	seen = { }
	bonds = [ seen.setdefault(b, b) for a in linker for b in a.bonds if b not in seen ]

	## Ligand atoms are new entities now
	ligand = cbase[0].residue.atoms
	static = cbase[0]
	## Set rotations
	for i, degrees in enumerate(ind[2]):
		try: 
			hyde5.createRotation(bonds[i], static)
			hyde5.rotate(bonds[i], degrees, absolute=True)
		except IndexError:
			break
	hbonds = hyde5.countHBonds(model, cache=False)
	clashes, num_of_clashes = hyde5.countClashes(atoms=ligand)
	hyde5.clearRotation(allbonds=True)

	## Set rotamers
	for i, rotId in enumerate(ind[3]):
		rotId = int(rotId)
		Rotamers.useRotamer(residues[i],[rotamers[i][rotId]])

	# Rotamer atoms are new entities now!
	#ligand_atoms = [ a for r in ligand for a in r.atoms ]
	all_atoms = [ a for m in chimera.openModels.list() for a in m.atoms]
	res_atoms = [ a for r in residues for a in r.atoms ]
	not_ligand = list(set(all_atoms) - set(ligand) - set(res_atoms))
	
	clashes_r, num_of_clashes_r = hyde5.countClashes(atoms=res_atoms,
		test=not_ligand)

	
	return len(hbonds), num_of_clashes, num_of_clashes_r

def hetCxOnePoint(ind1, ind2):

	for i, row in enumerate(ind1):
		if not i: #ignore ligand, fragment building
			continue
		size = min(len(ind1[i]), len(ind2[i]))
		if size > 1:
			cxpoint = random.randint(1, size - 1)
			ind1[i][cxpoint:], ind2[i][cxpoint:] = ind2[i][cxpoint:], ind1[i][cxpoint:]

	return ind1, ind2

def hetMutation(ind, indpb):
	for i, row in enumerate(ind):
		if random.random() < ind:
			if not i: 
				continue
			elif i == 1:
				j = random.randint(0,len(row)-1)
				ind[i][j] = random.randint(0,3)
				
			elif i == 2:
				j = random.randint(0,len(row)-1)
				ind[i][j] = random.randint(0,360)
				
			elif i == 3:
				j = random.randint(0,len(row)-1)
				ind[i][j] = random.uniform(0,8)
							
	return ind,

##/ FUNCTIONS

# Initial checks
if "base" not in chimera.selection.savedSels:
	raise UserError("Define a selection with the static part, named 'base'.")


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

# Set Chimera params
chimera.selection.saveSel("initial")
model = chimera.openModels.list()
residues = chimera.selection.savedSels['rotamers'].residues()
res_atoms = [a for res in residues for a in res.atoms]
rotamers = [ Rotamers.getRotamers(res)[1] for res in residues ]

# building blocks
wd = os.path.dirname(os.path.realpath(sys.argv[0]))
linkers = lego.getMol2Files(wd +'/mol2/linkers/')
fragments = lego.getMol2Files(wd + '/mol2/fragments/')

###
# Genetic Algorithm
####
# define individual, population, etc
deap.creator.create("FitnessMax", deap.base.Fitness, weights=(1.0, -1.0, -1.0))
deap.creator.create("Individual", list, fitness=deap.creator.FitnessMax)

toolbox = deap.base.Toolbox()
toolbox.register("rand_angle", random.randint, 0, 359)
toolbox.register("rand_h", random.randint, 1, 3)
toolbox.register("rand_rotamer", random.randint, 0, 8)
toolbox.register("rand_linker", random.randint, 0, len(linkers)-1)
toolbox.register("rand_fragment", random.randint, 0, len(fragments)-1)

toolbox.register("molecule", deap.tools.initCycle, list,
	[toolbox.rand_linker, toolbox.rand_fragment], n=1)
toolbox.register("h", deap.tools.initRepeat, list,
	toolbox.rand_h, n=2)
toolbox.register("linker_rots", deap.tools.initRepeat, list,
	toolbox.rand_angle, 5)
# toolbox.register("fragment_rots", deap.tools.initRepeat, list,
# 	toolbox.rand_angle)
toolbox.register("rotamers", deap.tools.initRepeat, list,
	toolbox.rand_rotamer, n=len(residues))

toolbox.register("individual", deap.tools.initCycle, 
	deap.creator.Individual, [toolbox.molecule, toolbox.h, toolbox.linker_rots,
	toolbox.rotamers], n=1)
toolbox.register("population", deap.tools.initRepeat, 
	list, toolbox.individual)


toolbox.register("evaluate", evalCoord)
toolbox.register("mate", hetCxOnePoint)
toolbox.register("mutate", hetMutation, indpb=0.05)
toolbox.register("select", deap.tools.selNSGA2)


def main():
	pop = toolbox.population(n=args.pop-1)
	hof = deap.tools.HallOfFame(1)
	stats = deap.tools.Statistics(lambda ind: ind.fitness.values)
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
    chimera.selection.setCurrent(chimera.selection.savedSels["initial"])
    evalCoord(hof[0])
    print("Best individual is: %s\nwith fitness: %s" % (hof[0], hof[0].fitness))
