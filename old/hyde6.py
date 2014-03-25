###########################################################
# CHIMERA SCRIPT
# Rotate selected bonds, find clashes and H bonds
# By Jaime RG <jaime.rogue@gmail.com>, in UAB
###########################################################

## VERSION 6
# Implement genetic algorithm

# TODO
# - Consider Von Mises Distribution for random angle
# - Calculate h bonds only if clashes < threshold?

import chimera, Rotamers, random, numpy, deap, argparse
from chimera import UserError
import hyde5, lego # the scripts!

# Initial checks
if "ligand" not in chimera.selection.savedSels:
	raise UserError("Define a selection named 'ligand' that\
	 contains the part of the molecule which will be moved")
elif len(chimera.selection.currentAtoms()) != 1:
	raise UserError("Select and anchor atom")
elif not len(chimera.selection.currentBonds()):
	raise UserError("Select some bonds to rotate")

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
bonds = chimera.selection.currentBonds()
anchor = chimera.selection.currentAtoms()[0]
model = chimera.selection.currentMolecules()
residues = chimera.selection.savedSels['rotamers'].residues()
res_atoms = [a for res in residues for a in res.atoms]
rotamers = [ Rotamers.getRotamers(res)[1] for res in residues ]
ligand = chimera.selection.savedSels['ligand'].atoms()
# TODO It will change with modifications!
#ligand_atoms = [ a for r in ligand for a in r.atoms ]
all_atoms = [ a for m in chimera.openModels.list() for a in m.atoms]
not_ligand = list(set(all_atoms) - set(ligand) - set(res_atoms))


for bond in bonds:
	hyde5.createRotation(bond, anchor)

def evalCoord(individual):
	
	for i, degrees in enumerate(individual[:len(bonds)]):
		hyde5.rotate(bonds[i], degrees, absolute=True)
	hbonds = hyde5.countHBonds(model)
	clashes, num_of_clashes = hyde5.countClashes(atoms=ligand)

	for i, rotId in enumerate(individual[len(bonds):]):
		rotId = int(rotId)
		Rotamers.useRotamer(residues[i],[rotamers[i][rotId]])

	res_atoms = [ a for r in residues for a in r.atoms ]
	clashes_r, num_of_clashes_r = hyde5.countClashes(atoms=res_atoms,
		test=not_ligand)

	
	return len(hbonds), num_of_clashes, num_of_clashes_r

def hetCxOnePoint(ind1, ind2, bound):

	ind1_0, ind1_1 = ind1[:bound], ind1[bound:]
	ind2_0, ind2_1 = ind2[:bound], ind2[bound:]
	size_0 = min(len(ind1_0), len(ind2_0))
	size_1 = min(len(ind1_1), len(ind2_1))

	if size_0 > 1:
		cxpoint = random.randint(1, size_0 - 1)
		ind1_0[cxpoint:], ind2_0[cxpoint:] = ind2_0[cxpoint:], ind1_0[cxpoint:]
	if size_1 > 1:
		cxpoint = random.randint(1, size_1 - 1)
		ind1_1[cxpoint:], ind2_1[cxpoint:] = ind2_1[cxpoint:], ind1_1[cxpoint:]

	ind1[:] = ind1_0 + ind1_1
	ind2[:] = ind2_0 + ind2_1

	return ind1, ind2


###
# Genetic Algorithm
####
# define individual, population, etc
deap.creator.create("FitnessMax", deap.base.Fitness, weights=(1.0, -1.0, -1.0))
deap.creator.create("Individual", list, fitness=deap.creator.FitnessMax)

toolbox = deap.base.Toolbox()
toolbox.register("rand_angle", random.randint, 0, 359)
toolbox.register("rand_rotamer", random.randint, 0, 8)
toolbox.register("individual", deap.tools.initCycle, 
	deap.creator.Individual, [toolbox.rand_angle] * len(bonds) + 
		[toolbox.rand_rotamer] * len(residues), n=1)
toolbox.register("population", deap.tools.initRepeat, 
	list, toolbox.individual)

LOW_BOUND = [0]*(len(bonds)+len(residues))
UP_BOUND = [360]*len(bonds) + [8]*len(residues)

toolbox.register("evaluate", evalCoord)
toolbox.register("mate", hetCxOnePoint, bound=len(bonds))
toolbox.register("mutate", deap.tools.mutPolynomialBounded,
	eta = 20.0, low = LOW_BOUND, up = UP_BOUND, indpb=0.05)
toolbox.register("select", deap.tools.selNSGA2)
adam = deap.creator.Individual([0] * len(bonds))

def main():
	pop = toolbox.population(n=args.pop-1) + [adam]
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
