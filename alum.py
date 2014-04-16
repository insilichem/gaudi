# aluminium optimization
# Rotamer optimization
import chimera, deap, Rotamers, mof3d
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
parser.add_argument('-t', '--threshold',
					required=False,
					default=2.0, 
					dest="threshold",
					metavar="<distance in angstroms | covalent>",
					help="Target distance for rotamers. If `covalent`, distance will be inferred from atom types." )
parser.add_argument('-w', '--wall',
					required=False,
					default=False, 
					dest="wall",
					action="store_true",
					help="Wether to allow distances < threshold or not" )
args = parser.parse_args()

#random.seed(64)

mol = chimera.openModels.list()[0]
residues = chimera.selection.savedSels['rotamers'].residues()
rotamers = [ Rotamers.getRotamers(r)[1] for r in residues ]
r_atoms = [ a for r in residues for a in r.atoms ]
al = chimera.selection.savedSels['al'].atoms()[0]

def evalCoord(ind):
	# rotamers
	for res, rot, i in zip(residues, rotamers, ind):
		try:
			Rotamers.useRotamer(res, [rot[i]])
		except IndexError:
			Rotamers.useRotamer(res, [rot[-1]])
	
	# distance
	avg_dist = mof3d.score.target.distance(residues, "O", al, args.threshold, \
		wall=args.wall, avg=True)
	# clashes
	r_atoms = [ a for r in residues for a in r.atoms ]
	clashes, num_of_clashes = mof3d.score.chem.clashes(atoms=r_atoms, 
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
	hof = deap.tools.HallOfFame(10)
	stats = deap.tools.Statistics(lambda ind: ind.fitness.values[0])
	stats.register("avg", numpy.mean, axis=0)
	stats.register("min", numpy.min, axis=0)
	stats.register("max", numpy.max, axis=0)
	pop, log = deap.algorithms.eaMuPlusLambda(pop, toolbox, 
		mu = int(args.pop/2), lambda_= int(args.pop/2), cxpb=0.5, 
		mutpb=0.2, ngen=args.ngen, stats=stats, halloffame=hof)
	return pop, log, hof

if __name__ == "__main__":

	pop, log, hof = main()
	evalCoord(hof[0])
	print("Best individual is: %s\nwith fitness: %s" % (hof[0], hof[0].fitness))
	print("More possible solutions to assess: ")
	for h in hof[1:]:
		print h, h.fitness