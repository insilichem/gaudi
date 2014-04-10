# Assess known crystals
# Rotamers only
import chimera
import deap
from deap import creator, algorithms, tools, base
import random, argparse, numpy, math
import hyde5
import ChemGroup as cg 
reload(hyde5)
def parseClashes(clashes):

	positive, negative = [], []
	for a1, c in clashes.items():
		for a2, dist in c.items():
			if dist <= 0.4 and a1.residue != a2.residue:
				if a1 in AROMATIC and a2 in AROMATIC:
					positive.append([a1, a2, dist, lj(a1, a2)])
				elif a1 in AROMATIC and a2 in ALIPHATIC:
					positive.append([a1, a2, dist, lj(a1, a2)])
				elif a1 in ALIPHATIC and a2 in AROMATIC:
					positive.append([a1, a2, dist, lj(a1, a2)])
			elif dist > 0.4:
				negative.append([a1, a2, dist, vdw_overlap(a1,a2)])

	return positive, negative

def lj(a1, a2):
	dist = a1.xformCoord().distance(a2.xformCoord())
	zero = 0.98*(a1.radius + a2.radius)
	x = zero/dist
	return (x**12 - 2*x**6)

def vdw_overlap(a1, a2):
	# Adapted from Eran Eyal, Comput Chem 25: 712-724, 2004
	d = a1.xformCoord().distance(a2.xformCoord())
	h_a, h_b = 0, 0
	if d < a1.radius+a2.radius:
		h_a = (a2.radius**2 - (d- a1.radius)**2)/(2*d)
		h_b = (a1.radius**2 - (d- a2.radius)**2)/(2*d)
	v = (math.pi/3) * (h_a**2) * (3*a1.radius - h_a) + \
		(math.pi/3) * (h_b**2) * (3*a2.radius - h_b)
	#maxv = (4/3.)*math.pi*(min(a1.radius, a2.radius)**3)
	return v

def evalCoord(ind, final=False):

	## 1 - Set rotations
	for i, br in enumerate(bondrots1+bondrots2):
		br.adjustAngle(ind[i] - br.angle, br.myanchor)
	
	## 2 - Score
	# TODO: Restrict donor and acceptors to smaller selection

	# hbonds_atoms = [ a for a in ligand_atoms if a.name not in ("C", "CA", "N", "O") ]
	# hbonds = hyde5.countHBonds([mol], sel=hbonds_atoms, cache=True)
	contacts, num_of_contacts =  hyde5.countClashes(
								atoms=ligand_atoms, test=mol.atoms, 
								intraRes=True, clashThreshold=-1.4, 
								hbondAllowance=0.0)

	positive_vdw, negative_vdw = parseClashes(contacts)


	if final:
		pbpos = chimera.misc.getPseudoBondGroup("hydrophobic interactions")
		pbneg = chimera.misc.getPseudoBondGroup("clashes")
		if positive_vdw: 
			max_pos, min_pos = max(abs(_[3]) for _ in positive_vdw), min(abs(_[3]) for _ in positive_vdw)
			out = open('/home/jr/x/positive_clashes.txt', 'w+')
			for p in positive_vdw:
				np = pbpos.newPseudoBond(p[0], p[1])
				intensity = (max_pos - abs(p[3]))/(max_pos - min_pos)
				opacity = 1 - 0.7*intensity
				np.color = chimera.MaterialColor(intensity,1,0,opacity)
				out.write("{0}\t{1}\n".format(p[2],p[3]))
			out.close()
		if negative_vdw:
			out = open('/home/jr/x/negative_clashes.txt', 'w+')
			max_neg, min_neg = max(_[2] for _ in negative_vdw), min(_[2] for _ in negative_vdw)
			for p in negative_vdw:
				np = pbneg.newPseudoBond(p[0], p[1])
				opacity = 1 - 0.7*abs(max_neg - p[2])/(max_neg - min_neg)
				np.color = chimera.MaterialColor(1,0,0,opacity)
				out.write("{0}\t{1}\n".format(p[2],p[3]))
			out.close()	
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

aliph = ['C3', [[cg.C, [cg.C, [cg.C , [cg.R, cg.R, cg.R, cg.R]], cg.R, cg.R] ], cg.R, cg.R, cg.R] ], [ 1, 1, 1, 1, 1, 0, 0] 
AROMATIC = set( a for g in cg.findGroup("aromatic ring", [mol]) for a in g )
ALIPHATIC = set( a for g in cg.findGroup(aliph, [mol]) for a in g if a not in AROMATIC )

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
