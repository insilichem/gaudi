# Assess known crystals
# Rotamers only
import chimera
import deap
from deap import creator, algorithms, tools, base
import random, argparse, numpy
import gaudi

def evalCoord(ind, final=False):

	## 1 - Set rotations
	for i, br in enumerate(bondrots1+bondrots2):
		br.adjustAngle(ind[i] - br.angle, br.myanchor)
	
	## 2 - Score
	# TODO: Restrict donor and acceptors to smaller selection

	# hbonds_atoms = [ a for a in ligand_atoms if a.name not in ("C", "CA", "N", "O") ]
	# hbonds = hyde5.countHBonds([mol], sel=hbonds_atoms, cache=True)
	contacts, num_of_contacts, positive_vdw, negative_vdw = \
		gaudi.score.chem.clashes(
								atoms=ligand_atoms, test=mol.atoms, 
								intraRes=True, clashThreshold=-1.4, 
								hbondAllowance=0.0, parse=True)
	if final:
		if positive_vdw: 
			gaudi.score.chem.draw_clashes(positive_vdw, startCol="FFFF00", endCol="00FF00",
				key=3, name="Hydrophobic interactions")
		if negative_vdw:
			gaudi.score.chem.draw_clashes(negative_vdw, startCol="FF0000", endCol="FF0000",
				key=3, name="Bad clashes")
	return  sum(abs(a[3]) for a in negative_vdw)/2, sum(1-a[3] for a in positive_vdw)/2, \
			0 #len(hbonds)

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

anchor1 = chimera.selection.savedSels['anchor1'].atoms()[0]
bondrots1 = []
for b in chimera.selection.savedSels['bonds1'].bonds():
	br = chimera.BondRot(b)
	br.myanchor = gaudi.utils.box.find_nearest(anchor1, b.atoms)
	bondrots1.append(br)

anchor2 = chimera.selection.savedSels['anchor2'].atoms()[0]
bondrots2 = []
for b in chimera.selection.savedSels['bonds2'].bonds():
	br = chimera.BondRot(b)
	br.myanchor = gaudi.utils.box.find_nearest(anchor2, b.atoms)
	bondrots2.append(br)

mol = chimera.openModels.list()[0]
ligands = chimera.selection.savedSels['ligands'].residues()
ligand_atoms = [ a for r in ligands for a in r.atoms ]


###
# Genetic Algorithm
# define individual, population, etc
deap.creator.create("FitnessMax", deap.base.Fitness, weights=(-1.0, 1.0, 1.0))
deap.creator.create("Individual", list, fitness=deap.creator.FitnessMax)

toolbox = deap.base.Toolbox()
toolbox.register("rand_angle", random.uniform, 0, 360)
toolbox.register("individual", deap.tools.initRepeat, deap.creator.Individual,
	toolbox.rand_angle, n=len(bondrots1)+len(bondrots2))
toolbox.register("population", deap.tools.initRepeat, 
	list, toolbox.individual)

toolbox.register("evaluate", evalCoord)
toolbox.register("mate", deap.tools.cxSimulatedBinaryBounded,
	eta=10., low=0., up=360.)
toolbox.register("mutate", deap.tools.mutPolynomialBounded, 
	eta=1., low=0., up=360., indpb=0.05)
toolbox.register("select", deap.tools.selNSGA2)

def main():
	numpy.set_printoptions(precision=4)
	pop = toolbox.population(n=args.pop)
	hof = deap.tools.HallOfFame(1)
	stats = deap.tools.Statistics(lambda ind: ind.fitness.values)
	stats.register("avg", numpy.mean, axis=0)
	stats.register("min", numpy.min, axis=0)
	stats.register("max", numpy.max, axis=0)
	try:
		pop, log = deap.algorithms.eaMuPlusLambda(pop, toolbox, 
			mu = int(args.pop/2), lambda_= int(args.pop/2), cxpb=0.5, 
			mutpb=0.2, ngen=args.ngen, stats=stats, halloffame=hof)
	except KeyboardInterrupt:
		return [], [], hof
	return pop, log, hof

if __name__ == "__main__":
	pop, log, hof = main()
	evalCoord(hof[0], final=True)
	print("Best individual is: %s\nwith fitness: %s" % (hof[0], hof[0].fitness))
