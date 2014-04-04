# Assess known crystals
# Rotamers only
import chimera
import deap
from deap import creator, algorithms, tools, base
import random, argparse, numpy
import hyde5
import ChemGroup as cg 
reload(hyde5)
def parseClashes(clashes):
	aromatoms = set( a for g in cg.findGroup("aromatic ring", [mol]) for a in g )
	positive, negative = [], []
	for a1, c in clashes.items():
		for a2, dist in c.items():
			if a1 in aromatoms and a2 in aromatoms and dist<=0.4:
				positive.append([a1, a2, dist])
			elif dist >= 0.6:
				negative.append([a1, a2, dist])

	return positive, negative

def evalCoord(ind):

	## 1 - Set rotations
	for i, br in enumerate(bondrots1+bondrots2):
		br.adjustAngle(ind[i] - br.angle, br.myanchor)
	
	## 2 - Score
	# TODO: Restrict donor and acceptors to smaller selection

	hbonds_atoms = [ a for a in ligand_atoms if a.name not in ("C", "CA", "N", "O") ]
	hbonds = hyde5.countHBonds([mol], sel=hbonds_atoms, cache=True)
	contacts, num_of_contacts =  hyde5.countClashes(
								atoms=ligand_atoms, test=mol.atoms, 
								intraRes=True, clashThreshold=-0.4, 
								hbondAllowance=0.0)

	positive_contacts, negative_contacts = parseClashes(contacts)
	chimera.selection.setCurrent([ l[0] for l in positive_contacts if l[0] in ligand_atoms ])
	return len(negative_contacts)/2, len(positive_contacts)/2, len(hbonds)


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
	br.myanchor = hyde5.findNearest(anchor1, b.atoms)
	bondrots1.append(br)

anchor2 = chimera.selection.savedSels['anchor2'].atoms()[0]
bondrots2 = []
for b in chimera.selection.savedSels['bonds2'].bonds():
	br = chimera.BondRot(b)
	br.myanchor = hyde5.findNearest(anchor2, b.atoms)
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
