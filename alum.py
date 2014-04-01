# aluminium optimization
# Rotamer optimization
import chimera, deap, Rotamers, hyde5
from deap import creator, algorithms, tools, base
import random, argparse, numpy

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

#random.seed(64)

mol = chimera.openModels.list()[0]
residues = chimera.selection.savedSels['rotamers'].residues()
rotamers = [ Rotamers.getRotamers(r)[1] for r in residues ]
r_atoms = [ a for r in residues for a in r.atoms ]
al = chimera.selection.savedSels['al'].atoms()[0]

def evalCoord(ind):
	for i, r in enumerate(ind):
		try:
			Rotamers.useRotamer(residues[i], [rotamers[i][r]])
		except IndexError:
			Rotamers.useRotamer(residues[i], [rotamers[i][-1]])
	# distance
	r_atoms = [ a for r in residues for a in r.atoms ]
	
	distances = []
	for r in residues:
		oxygens = [ a for a in r.atoms if a.element.number == 8]
		d = [ ox.xformCoord().distance(al.xformCoord()) for ox in oxygens ]
		# if all([d_ > 2.0 for d_ in d]):
		# 	distances.append(min(d, key = lambda x: abs(x-2.0)))
		# else:
		# 	distances.append(10.0)
		distances.append(min(d, key = lambda x: abs(x-2.0)))
		
	avg_dist = numpy.mean(distances)
	clashes, num_of_clashes = hyde5.countClashes(atoms=r_atoms, 
		test=[ a for a in mol.atoms if a != al ])

	return avg_dist, num_of_clashes

# Genetic Algorithm
# define individual, population, etc
deap.creator.create("FitnessMin", deap.base.Fitness, weights=(-1.0, -1.0))
deap.creator.create("Individual", list, fitness=deap.creator.FitnessMin)

toolbox = deap.base.Toolbox()
toolbox.register("randrot", random.randint, 0, 8)
toolbox.register("individual", deap.tools.initRepeat, deap.creator.Individual,
	toolbox.randrot, n=len(residues))
toolbox.register("population", deap.tools.initRepeat, 
	list, toolbox.individual)

toolbox.register("evaluate", evalCoord)
toolbox.register("mate", deap.tools.cxTwoPoint)
toolbox.register("mutate", deap.tools.mutUniformInt, 
	low=0., up=8., indpb=0.05)
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