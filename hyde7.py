###########################################################
# CHIMERA SCRIPT
# Rotate selected bonds, find clashes and H bonds
# By Jaime RG <jaime.rogue@gmail.com>, in UAB
###########################################################

## VERSION 7
# Implement genetic algorithm

# TODO
# - Calculate h bonds only if clashes < threshold?


# Chimera
import chimera, Rotamers, SwapRes, SplitMolecule, ChemGroup as cg
from chimera import UserError
# Python
import random, numpy, deap, argparse, sys, os
from deap import creator, tools, base, algorithms
# Custom
import hyde5, lego
import fragment3 as frag
from itertools import product
### CUSTOM FUNCTIONS

def molLibrary(cbase, linkers, fragments, dihedral=None, alpha=None):
	
	explore = product(range(len(linkers)), range(len(fragments)))
	library = {}
	for x in explore:
		i, j = x
		new = SplitMolecule.split.molecule_from_atoms(mol, cbase.atoms)
		target = new.atoms[-1]
		linker = frag.insertMol(linkers[i], target=target, alpha=alpha, alpha2=120.0)
		linker_anchor = [ a for a in linker if a.anchor in (4,6,8)][0]
		frag.insertMol(fragments[j], target=linker_anchor, alpha2=114.125)
		
		fragment_anchor = [ a for a in linker_anchor.neighbors if a.element.number != 1
			and a not in linker ]
		target_neighbor = [ a for a in target.neighbors if a.element.number != 1
			and a not in linker ][0]
		bonds = getSequentialBonds([target]+linker[:]+fragment_anchor,target_neighbor)
		
		bondrots = []
		for b in bonds:
			br = chimera.BondRot(b)
			br.myanchor = hyde5.findNearest(new.atoms[0], b.atoms)
			bondrots.append(br)

		library[i,j] = [new, bondrots]
	return library

def getSequentialBonds(atoms,s):
	if s not in atoms: atoms.append(s)
	bonds = list(set([ b for a in atoms for b in a.bonds if set(b.atoms)<=set(atoms) ]))
	nbonds = []
	while bonds:
		b = bonds.pop(0)
		if s in b.atoms:
			nbonds.append(b)
			s = b.otherAtom(s)
		else: 
			bonds.append(b)
	return nbonds

def parseClashes(clashes):
	aliph = ['C3', [[cg.C, [cg.C, [cg.C , [cg.R, cg.R, cg.R, cg.R]], cg.R, cg.R] ], cg.R, cg.R, cg.R] ], [ 1, 1, 1, 1, 1, 0, 0] 
	AROMATIC = set( a for g in cg.findGroup("aromatic ring", [mol]) for a in g )
	ALIPHATIC = set( a for g in cg.findGroup(aliph, [mol]) for a in g if a not in AROMATIC )

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
	v = (numpy.pi/3) * (h_a**2) * (3*a1.radius - h_a) + \
		(numpy.pi/3) * (h_b**2) * (3*a2.radius - h_b)
	#maxv = (4/3.)*numpy.pi*(min(a1.radius, a2.radius)**3)
	return v

def evalCoord(ind, close=True, hidden=False):
	## 1 - Choose ligand from pre-built mol library
	ligand, bondrots = mol_library[ind['molecule'][0],ind['molecule'][1]]
	chimera.openModels.add([ligand], shareXform=True, hidden=hidden)
	# Chimera converts metal bonds to pseudoBonds all the time
	pbgroup = chimera.misc.getPseudoBondGroup(
					"coordination complexes of %s (%s)" % 
					(ligand.name, ligand), associateWith=[ligand])
	if pbgroup.pseudoBonds:
		for pb in pbgroup.pseudoBonds:
			chimera.molEdit.addBond(pb.atoms[0],pb.atoms[1])
		pbm = ligand.pseudoBondMgr()
		pbm.deletePseudoBondGroup(pbgroup)

	## 2 - Set rotations
	# Direct access to BondRot, instead of BondRotMgr
	for i, br in enumerate(bondrots):
		br.adjustAngle(ind['linker_rots'][i] - br.angle, br.myanchor)

	## 3 - Set mutamers/rotamers
	for i, aa in enumerate(ind['mutamers']):
		try: 
			rotamers = Rotamers.getRotamers(residues[i], resType=aminoacids[aa])[1] 
			rotId = ind['rotamers'][i]
			Rotamers.useRotamer(residues[i],[rotamers[rotId]])
		except Rotamers.NoResidueRotamersError: # ALA, GLY...
			SwapRes.swap(residues[i], aminoacids[aa], bfactor=None)
		except IndexError:
			Rotamers.useRotamer(residues[i],[rotamers[-1]])

	## 4 - Score
	# Inserted atoms become new entities
	res_atoms = [ a for r in residues for a in r.atoms ]
	models = chimera.openModels.list(all=True,modelTypes=[chimera.Molecule])
	ligand_env.clear()
	ligand_env.add(ligand)
	ligand_env.merge(chimera.selection.REPLACE, 
					chimera.specifier.zone(ligand_env, 'atom', None, 15.0, models))

	hbonds = hyde5.countHBonds(
				models, cache=False, test=ligand_env.atoms(),
				sel=[ a for a in ligand.atoms if a not in ("C", "CA", "N", "O") ])
	# TODO: Restrict test to smaller selection
	contacts, num_of_contacts =  hyde5.countClashes(
									atoms=ligand.atoms, 
									test=ligand_env.atoms(), 
									intraRes=True, clashThreshold=-0.4, 
									hbondAllowance=0.0)
	chimera.selection.setCurrent(ligand_env.atoms())
	clashes_r, num_of_clashes_r = hyde5.countClashes(
									atoms=res_atoms,
									test=mol.atoms)
	
	positive_vdw, negative_vdw = parseClashes(contacts)

	if not close:
		pbpos = chimera.misc.getPseudoBondGroup("hydrophobic interactions")
		pbneg = chimera.misc.getPseudoBondGroup("clashes")
		if positive_vdw: 
			max_pos= max(abs(_[3]) for _ in positive_vdw)
			for p in positive_vdw:
				np = pbpos.newPseudoBond(p[0], p[1])
				intensity = (max_pos - abs(p[3]))/(max_pos)
				opacity = 1 - 0.7*intensity
				np.color = chimera.MaterialColor(intensity,1,0,opacity)
		if negative_vdw:
			max_neg, min_neg = max(_[2] for _ in negative_vdw), min(_[2] for _ in negative_vdw)
			for p in negative_vdw:
				opacity = 1
				np = pbneg.newPseudoBond(p[0], p[1])
				if max_neg != min_neg:
					opacity = 1 - 0.7*abs(max_neg - p[2])/(max_neg - min_neg)
				np.color = chimera.MaterialColor(1,0,0,opacity)
	else:
		chimera.openModels.remove([ligand])

	return  len(hbonds), sum(abs(a[3]) for a in negative_vdw)/2, \
			sum(1-a[3] for a in positive_vdw)/2, num_of_clashes_r

def hetCxOnePoint(ind1, ind2):
	
	for key in ind1:
		if key == 'molecule': 
			continue
		elif key == 'linker_rots':
			ind1[key][:], ind2[key][:] = deap.tools.cxSimulatedBinaryBounded(
				ind1[key], ind2[key], eta=10., low=0., up=360.)
		elif key == 'mutamers':
			ind1[key], ind2[key] = deap.tools.cxTwoPoint(ind1[key], ind2[key])
		elif key == 'rotamers':
			ind1[key], ind2[key] = deap.tools.cxTwoPoint(ind1[key], ind2[key])
	return ind1, ind2

def hetMutation(ind, indpb):
	
	for key in ind:
		if key == 'molecule': 
			continue
		elif key == 'linker_rots':
			ind[key] = deap.tools.mutPolynomialBounded(ind[key], 
				eta=10., low=0., up=360., indpb=indpb)[0]
		elif key == 'mutamers':
			ind[key] = deap.tools.mutUniformInt(ind[key], 
				low=0, up=len(residues)-1, indpb=indpb)[0]
		elif key == 'rotamers':
			ind[key] = deap.tools.mutUniformInt(ind[key], 
				low=0, up=8, indpb=indpb)[0]
	return ind,

# Couple of constants
aminoacids = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLU', 'GLN',
			  'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE',
			  'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL' ]
##/ FUNCTIONS

# Initial checks
required_sels = 'base', 'anchor', 'mutable'
if not set(required_sels) <= set(chimera.selection.savedSels):
	raise UserError("""You have to define three selections with 'namesel'
---
- 'base': Single atom. Terminal end of the static part of the ligand
- 'anchor': Single atom. Together with 'base', it will define the static part. It's also the insertion point for the linker.
- 'mutable': Residues to be swapped and mutated.
""")

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
parser.add_argument('-e', '--eval',
					required=False,
					type=str,
					default='', 
					dest="eval",
					metavar="<individual>",
					help="Paste individual data to evaluate" )
parser.add_argument('-pf', '--pareto',
					required=False,
					default=False, 
					dest="pareto",
					action="store_true",
					help="Wether to report Pareto Front rank or not" )
args = parser.parse_args()
wd = os.path.dirname(os.path.realpath(sys.argv[0]))

# Get Chimera params
residues = [ r for r in sorted(chimera.selection.savedSels['mutable'].residues(), 
	key= lambda r: (r.id.chainId, r.id.position)) ]
base_at = chimera.selection.savedSels['base'].atoms()[0]
mol = base_at.molecule
ligand = base_at.residue
anchor = chimera.selection.savedSels['anchor'].atoms()[0]
cbase = [base_at] + list(hyde5.atomsBetween(base_at, anchor)) + [anchor]
ligand_env = chimera.selection.ItemizedSelection()
#dihedral
d3 = anchor
d4 = [a for a in anchor.neighbors if a not in cbase][0]
d2 = [a for a in d3.neighbors if a in cbase][0]
d1 = [a for a in d2.neighbors if a != anchor ]
if len(d1)>1: 
	d1 = [a for a in d1 if a.numBonds>1]
d1 = d1[0]
dihedral = chimera.dihedral(d1.coord(),d2.coord(),d3.coord(),d4.coord())
#angle
alpha = chimera.angle(d2.coord(),d3.coord(),d4.coord())

# Save cbase in ghost molecule
cbasecopy = SplitMolecule.split.molecule_from_atoms(mol, cbase)
[mol.deleteAtom(a) for a in ligand.atoms]

# Get building blocks
linkers = sorted(lego.getMol2Files(wd +'/mol2/linkers/'))
fragments = sorted(lego.getMol2Files(wd + '/mol2/fragments/'))

# Build library
mol_library = molLibrary(cbasecopy, linkers, fragments, dihedral=dihedral, alpha=alpha)

###
# Genetic Algorithm
# define individual, population, etc
deap.creator.create("FitnessMax", deap.base.Fitness, weights=(1.0, -1.0, 1.0, -1.0))
deap.creator.create("Individual", dict, fitness=deap.creator.FitnessMax)

# Operators
toolbox = deap.base.Toolbox()
toolbox.register("rand_angle", random.uniform, 0, 360)
toolbox.register("rand_aa", random.randint, 0, len(aminoacids)-1)
toolbox.register("rand_rotamer", random.randint, 0, 8)
toolbox.register("rand_linker", random.randint, 0, len(linkers)-1)
toolbox.register("rand_fragment", random.randint, 0, len(fragments)-1)

# Genes
toolbox.register("molecule", deap.tools.initCycle, list,
	[toolbox.rand_linker, toolbox.rand_fragment], n=1)
toolbox.register("linker_rots", deap.tools.initRepeat, list,
	toolbox.rand_angle, n=8)
# toolbox.register("fragment_rots", deap.tools.initRepeat, list,
# 	toolbox.rand_angle)
toolbox.register("mutamers", deap.tools.initRepeat, list,
	toolbox.rand_aa, n=len(residues))
toolbox.register("rotamers", deap.tools.initRepeat, list,
	toolbox.rand_rotamer, n=len(residues))
toolbox.register("toDict", 
	(lambda ind, *fn: ind((f.__name__, f()) for f in fn)))

# Individual and population
toolbox.register("individual", toolbox.toDict, deap.creator.Individual, 
	toolbox.molecule, toolbox.linker_rots, toolbox.mutamers, toolbox.rotamers)
toolbox.register("population", deap.tools.initRepeat, list, toolbox.individual)

# Aliases for algorithm
toolbox.register("evaluate", evalCoord)
toolbox.register("mate", hetCxOnePoint)
toolbox.register("mutate", hetMutation, indpb=0.05)
toolbox.register("select", deap.tools.selNSGA2)

def main():
	pop = toolbox.population(n=args.pop)
	if args.pareto:
		hof = deap.tools.ParetoFront()
	else: 
		hof = deap.tools.HallOfFame(20)
	stats = deap.tools.Statistics(lambda ind: ind.fitness.values)
	stats.register("avg", numpy.mean, axis=0)
	stats.register("min", numpy.min, axis=0)
	stats.register("max", numpy.max, axis=0)
	pop, log = deap.algorithms.eaMuPlusLambda(pop, toolbox, 
		mu = int(args.pop/2), lambda_= int(args.pop/2), cxpb=0.5, 
		mutpb=0.2, ngen=args.ngen, stats=stats, halloffame=hof)
	return pop, log, hof

if __name__ == "__main__":
	numpy.set_printoptions(precision=2)
	if args.eval:
		print "Fitness: ", evalCoord(eval(args.eval), close=False)
	else:
		pop, log, hof = main()
		evalCoord(hof[0], close=False)
		print("Best individual is: %s\nwith fitness: %s" % (hof[0], hof[0].fitness))
		print("More possible solutions to assess:")
		for h in hof[1:11]:
			print h, h.fitness
