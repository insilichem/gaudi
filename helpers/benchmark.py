import deap, gaudi, os, glob, numpy, sys, chimera, random
from deap import creator, tools, base, algorithms
evalcount = 0
def evaluate(ind):
	savedir = ''.join(cfg.default.savepath, id(ind), '/')
	os.mkdir(savedir)
	weights = [ w*r for (w,r) in zip(cfg_weights, ind) ]
	print '\nBenchmarking weights', weights
	chimera.runCommand('runscript /home/jr/x/gaudi/base.py ' + \
						sys.argv[1] + ' ' + ' '.join(map(str,weights)))
	assess, = chimera.openModels.open(cfg.ligand.assess, shareXform=True)
	results = glob.glob(cfg.default.savepath+"*.mol2")
	rmsd, fitness = [], []
	for result in results:
		r, = chimera.openModels.open(result)
		fit = r.mol2comments[r.mol2comments.index('#>>GAUDI.score')+2].split()[1:]
		fitness.append(sum(w*f for (w,f) in zip(weights, fit)))
		rmsd.append(gaudi.utils.box.rmsd(r, assess))
		chimera.openModels.remove([r])
		os.rename(result, savedir+result)
	chimera.runCommand('close all')

	sort_by_rmsd = sorted(results, key=rmsd)
	sort_by_fitness = sorted(results, key=fitness)

	return numpy.mean(x==y for (x,y) in zip(sort_by_rmsd, sort_by_fitness)),

def main():
	pop = toolbox.population(n=int(sys.argv[2]))
	hof = deap.tools.ParetoFront()
	stats = deap.tools.Statistics(lambda ind: ind.fitness.values)
	stats.register("avg", numpy.mean, axis=0)
	stats.register("min", numpy.min, axis=0)
	stats.register("max", numpy.max, axis=0)
	pop, log = deap.algorithms.eaMuPlusLambda(pop, toolbox, 
		mu=int(0.75*int(sys.argv[2])), lambda_=int(0.75*int(sys.argv[2])), 
		cxpb=0.5, mutpb=0.25, 
		ngen=int(sys.argv[3]), stats=stats, halloffame=hof)
	return pop, log, hof


cfg = gaudi.utils.parse.Settings(sys.argv[1])
cfg_weights = cfg.weights()

toolbox = deap.base.Toolbox()
deap.creator.create("FitnessMin", deap.base.Fitness, weights=(-1.0,))
deap.creator.create("Individual", list, fitness=deap.creator.FitnessMin)

# Individual and population

toolbox.register("individual", deap.tools.initRepeat, deap.creator.Individual,
	random.random, n=len(cfg_weights))
toolbox.register("population", deap.tools.initRepeat, list, toolbox.individual)

# Aliases for algorithm
toolbox.register("evaluate", evaluate)
toolbox.register("mate", deap.tools.cxSimulatedBinaryBounded,
	eta=10., low=0, up=1)
toolbox.register("mutate", deap.tools.mutPolynomialBounded, 
	eta=10., low=0, up=1, indpb=0.05)
toolbox.register("select", deap.tools.selNSGA2)

if __name__ == "__main__":	
	pop, log, hof = main()
	print "\n---------------\n"
	for h in hof:
		print id(h), h, '\t', h.fitness

