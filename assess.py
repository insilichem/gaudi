# Assess known crystals
# Rotamers only
import chimera
import deap
from deap import creator, algorithms, tools, base
import random, argparse, numpy
import hyde5

def evalCoord(ind):

	## 1 - Set rotations
	for i, br in enumerate(bondrots):
		#print i, br.bond, br.myanchor
		chimera.selection.setCurrent(br.myanchor)
		br.adjustAngle(ind[i] - br.angle, br.myanchor)

	## 2 - Score
	# TODO: Restrict donor and acceptors to smaller selection
	model = chimera.openModels.list()
	hbonds = hyde5.countHBonds(model, sel=ligand.atoms, cache=True)
	clashes, num_of_clashes = hyde5.countClashes(atoms=ligand.atoms, 
		test=mol.atoms)

	return num_of_clashes, len(hbonds)


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

# Remove after testing
random.seed(64)

anchor = chimera.selection.currentAtoms()[0]
bondrots = []
for b in chimera.selection.currentBonds():
	br = chimera.BondRot(b)
	br.myanchor = hyde5.findNearest(anchor, b.atoms)
	bondrots.append(br)

mol = chimera.selection.currentMolecules()[0]
ligand = chimera.selection.currentResidues()[0]

###
# Genetic Algorithm
# define individual, population, etc
deap.creator.create("FitnessMax", deap.base.Fitness, weights=(-1.0, 1.0))
deap.creator.create("Individual", list, fitness=deap.creator.FitnessMax)

toolbox = deap.base.Toolbox()
toolbox.register("rand_angle", random.uniform, 0, 360)
toolbox.register("individual", deap.tools.initRepeat, deap.creator.Individual,
	toolbox.rand_angle, n=len(bondrots))
toolbox.register("population", deap.tools.initRepeat, 
	list, toolbox.individual)

toolbox.register("evaluate", evalCoord)
toolbox.register("mate", deap.tools.cxSimulatedBinaryBounded,
	eta=10., low=0., up=360.)
toolbox.register("mutate", deap.tools.mutPolynomialBounded, 
	eta=10., low=0., up=360., indpb=0.05)
toolbox.register("select", deap.tools.selNSGA2)

def main():
	pop = toolbox.population(n=args.pop)
	hof = deap.tools.HallOfFame(1)
	stats = deap.tools.Statistics(lambda ind: ind.fitness.values[0])
	stats.register("avg", numpy.mean)
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
