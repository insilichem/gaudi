# aluminium optimization
# Rotamer optimization
import chimera, deap, Rotamers, gaudi, Matrix as M
from deap import creator, algorithms, tools, base
import random, argparse, numpy, math

parser = argparse.ArgumentParser()
parser.add_argument('-p', '--population',
					required=False,
					type=int,
					default=300, 
					dest="pop",
					metavar="<number of individuals>",
					help="Number of individuals in each population" )
parser.add_argument('-g', '--generations',
					required=False,
					type=int,
					default=15, 
					dest="ngen",
					metavar="<number of generations>",
					help="Number of generations to calculate" )
parser.add_argument('-t', '--threshold',
					required=False,
					type=float,
					default=2.2, 
					dest="threshold",
					metavar="<distance in angstroms | `covalent`>",
					help="Target distance for rotamers. If `covalent`, distance will be inferred from atom types." )
parser.add_argument('-u', '--threshold2',
					required=False,
					default=1.8, 
					dest="threshold2",
					type=float,
					help="Wether to allow distances < threshold or not" )
parser.add_argument('-r', '--radius',
					required=False,
					default=2.5, 
					type=float,
					dest="radius",
					metavar="<angstroms>",
					help="Radius of the sphere where the metal will be able to move"  )
args = parser.parse_args()

################################

## Find aluminium
mol = chimera.openModels.list()[0]
alum = [a for r in mol.residues if r.type=='AL3' for a in r.atoms]
alcrd, alelem = alum[0].coord(), alum[0].element
# New movable aluminium
alres = gaudi.molecule._dummy_res('alum')
al = chimera.molEdit.addAtom('ALX', alelem, alres, alcrd, bondedTo=None, serialNumber=-1)
chimera.selection.setCurrent(al)
zone = chimera.specifier.evalSpec('sel zr < 7')
# Original aluminium ions are not needed
[mol.deleteResidue(r) for r in mol.residues if r.type=='AL3']
# Discover nearby residues
residues = [ r for r in zone.residues() if r.type in ('ASP', 'GLU') ]
rotamers = [ Rotamers.getRotamers(r)[1] for r in residues ]
r_atoms = [ a for r in residues for a in r.atoms ]

def evalCoord(ind, close=True):
	## Rebuild system
	# rotamers
	for res, rot, i in zip(residues, rotamers, ind['rotamers']):
		try:
			Rotamers.useRotamer(res, [rot[i]])
		except IndexError:
			Rotamers.useRotamer(res, [rot[-1]])
	# move metal
	al.molecule.openState.xform = \
			M.chimera_xform(M.multiply_matrices(ind['position'][0], ind['position'][2]))
	if not close:
		return mol, al.molecule
	
	## Test new conformation
	# distance
	zone = chimera.specifier.evalSpec('sel zr < 7')
	oxygens = [ (gaudi.score.target._distance(a, al),a)
			for a in zone.atoms() if a.element.number==8 and len(a.neighbors)==1 ]
	oxygens = [ ox for ox in oxygens if ox[0]>args.threshold2]
	oxygens.sort()

	# dihedral
	dihedrals = []
	for d,ox in oxygens[:3]:
		nearest_ox = ox
		nearest_ox_c = nearest_ox.neighbors[0]
		nearest_ox_cc = next(a for a in nearest_ox_c.neighbors if a.element.number!=8)
		dihedral = chimera.dihedral(*[a.molecule.openState.xform.apply(a.coord())
									for a in (al, nearest_ox, nearest_ox_c, nearest_ox_cc)])
		planarity = abs(math.sin(math.radians(dihedral)))
		dihedrals.append(planarity)

	# clashes
	r_atoms = [ a for r in residues for a in r.atoms ] + [al]
	clashes, num_of_clashes, pos, neg = gaudi.score.chem.clashes(atoms=r_atoms, 
		test=(a for a in mol.atoms), parse=True)

	score = [len([d for (d,o) in oxygens if d < args.threshold])]
	score.append(sum(abs(a[3]) for a in neg)/2)
	score.extend([abs(2.0-ox[0]) for ox in oxygens[:3]] + [0]*(3-len(oxygens[:3])))
	score.extend(dihedrals + [0]*(3-len(dihedrals)))
	
	return score


def het_crossover(ind1, ind2):
	ind1['rotamers'], ind2['rotamers'] = deap.tools.cxTwoPoint(ind1['rotamers'], ind2['rotamers'])
	return ind1, ind2

def het_mutation(ind, indpb):
	for key in ind:
		if key == 'rotamers':
			ind[key], = deap.tools.mutUniformInt(ind[key], 
				low=0, up=10, indpb=indpb)
		elif key == 'xform' and random.random() < indpb:
			# Careful! Mutation generates a whole NEW position (similar to eta ~= 0)
			ind['xform'] = gaudi.move.rand_xform(al[0].coord(), args.radius, rotate=False)

	return ind,

# Genetic Algorithm
# define individual, population, etc
deap.creator.create("FitnessMin", deap.base.Fitness, weights=(1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0))
deap.creator.create("Individual", dict, fitness=deap.creator.FitnessMin,
					fitness_names=[ 'nearby_ox', 'clashes', 'distance', 'distance', 'distance', 'planarity'])

toolbox = deap.base.Toolbox()
toolbox.register("toDict", 
	(lambda ind, *fn: ind((f.__name__, f()) for f in fn)))
#genes
toolbox.register("randrot", random.randint, 0, 10)
toolbox.register("rotamers", deap.tools.initRepeat, list, toolbox.randrot, n=len(residues))
toolbox.register("position", gaudi.move.rand_xform, alcrd, args.radius, rotate=False)
genes = [toolbox.rotamers, toolbox.position]
#ind&pop
toolbox.register("individual", toolbox.toDict, deap.creator.Individual, *genes)
toolbox.register("population", deap.tools.initRepeat, list, toolbox.individual)

toolbox.register("evaluate", evalCoord)
toolbox.register("mate", het_crossover)
toolbox.register("mutate", het_mutation, indpb=0.05)
toolbox.register("select", deap.tools.selNSGA2)

def main():
	pop = toolbox.population(n=args.pop)
	hof = deap.tools.ParetoFront()
	stats = deap.tools.Statistics(lambda ind: ind.fitness.values)
	stats.register("avg", numpy.mean, axis=0)
	stats.register("min", numpy.min, axis=0)
	stats.register("max", numpy.max, axis=0)
	pop, log = deap.algorithms.eaMuPlusLambda(pop, toolbox, 
		mu=int(args.pop/2), lambda_=int(args.pop/2), cxpb=0.5, 
		mutpb=0.2, ngen=args.ngen, stats=stats, halloffame=hof)
	return pop, log, hof

if __name__ == "__main__":
	pop, log, hof = main()
	rank = gaudi.utils.box.write_individuals(hof, mol.openedAs[0][:-4]+"/", 'solution', evalCoord, remove=False)
	out = open(mol.openedAs[0][:-4]+'/results.gaudi', 'w+')
	print >> out,  '#~ Generated by GAUDI\n'
	print >> out,  '>>GAUDI.results'
	print >> out,  '\n'.join(rank)
	out.close()