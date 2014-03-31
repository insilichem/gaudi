###########################################################
# CHIMERA SCRIPT
# Rotate selected bonds, find clashes and H bonds
# By Jaime RG <jaime.rogue@gmail.com>, in UAB
###########################################################

## VERSION 7
# Implement genetic algorithm

# TODO
# - Consider Von Mises Distribution for random angle
# - Calculate h bonds only if clashes < threshold?
# - Avoid ligand rebuilding every step! Use 'ghost' library instead
# - Better crossover and mutations functions

# Chimera
import chimera, Rotamers, SwapRes, SplitMolecule
from chimera import UserError
# Python
import random, numpy, deap, argparse, sys, os
from deap import creator, tools, base, algorithms
# Custom
import hyde5, lego
import fragment3 as frag
from itertools import product
reload(frag)
### CUSTOM FUNCTIONS

def molLibrary(cbase, linkers, fragments, link_end=3, dihedral=120.0, alpha=120.0):
	
	explore = product(range(len(linkers)), range(len(fragments)))
		#range(base_geom), range(link_end))
	
	h1, h2 = 1, 2 #this will be part of GA
	library = {}
	for x in explore:
		i, j = x
		new = SplitMolecule.split.molecule_from_atoms(mol, cbase.atoms)
		target = new.atoms[-1]
		linker = frag.insertMol(linkers[i], target=target, h=h1, 
			alpha=alpha, dihedral=dihedral)
		linker_anchor = [ a for a in linker if a.anchor in (4,6,8)][0]
		frag.insertMol(fragments[j], target=linker_anchor, h=h2, alpha=-120.0)
		
		fragment_anchor = [ a for a in linker_anchor.neighbors if a.element.number != 1
			and a not in linker ]
		# chimera.openModels.add([new], shareXform=True)
		# chimera.selection.setCurrent([target]+linker[:]+fragment_anchor)
		bonds = getSequentialBonds(linker[:]+fragment_anchor,target)
		
		bondrots = []
		for b in bonds:
			br = chimera.BondRot(b)
			br.myanchor = hyde5.findNearest(new.atoms[0], b.atoms)
			bondrots.append(br)

		# Include h1, h2 in the future
		# library[i,j,h1,h2] = [new, bonds]
		library[i,j] = [new, bonds, bondrots]
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

def evalCoord(ind, close=True):

	# # 1 - Build ligand (DEPRECATED; It built the ligand every time
	# lego.clearBase(cbase)
	# linker = frag.insertMol(linkers[ind['molecule'][0]], target=anchor, join=True, 
	# 	inplace=True, h=ind['h1'][0])
	# linker_anchor = [ a for a in linker if a.anchor in (4,6,8) ]
	# frag.insertMol(fragments[ind['molecule'][1]], target=linker_anchor[0], 
	# 	alpha=-120., h=ind['h2'][0])		
	# seen = { } # Used to remove duplicate entries in `bonds`
	# bonds = [ seen.setdefault(b, b) for a in linker for b in a.bonds if b not in seen ]

	## 1 - Choose ligand from pre-built mol library
	ligand, bonds, bondrots = mol_library[ind['molecule'][0],ind['molecule'][1]]#,ind['h1'][0],ind['h2'][0]]
	chimera.openModels.add([ligand], shareXform=True)
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
		#hyde5.bondrot(bond, anchor=ligand.atoms[0], delta=ind['linker_rots'][i])
		#br = bondrots[i]
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
	# TODO: Restrict donor and acceptors to smaller selection
	model = chimera.openModels.list()
	hbonds = hyde5.countHBonds(model, sel=ligand.atoms, cache=False)
	# TODO: Restrict test to smaller selection
	clashes, num_of_clashes = hyde5.countClashes(atoms=ligand.atoms, 
		test=mol.atoms + ligand.atoms) 
	clashes_r, num_of_clashes_r = hyde5.countClashes(atoms=res_atoms,
		test=mol.atoms)
	
	if close:
		chimera.openModels.remove([ligand])

	return len(hbonds), num_of_clashes, num_of_clashes_r

def hetCxOnePoint(ind1, ind2):

	for key in ind1:
		if key in ('molecule', 'h1', 'h2'): #ignore building blocks FOR NOW ;)
			continue
		size = min(len(ind1[key]), len(ind2[key]))
		if size > 1:
			cxpoint = random.randint(1, size - 1)
			ind1[key][cxpoint:], ind2[key][cxpoint:] = \
			ind2[key][cxpoint:], ind1[key][cxpoint:]

	return ind1, ind2

def hetMutation(ind, indpb):
	for key, row in ind.items():
		if random.random() < ind:
			j = random.randint(0,len(row)-1)
			if key == 'molecule': 
				continue
			elif key == 'h1':
				ind[key][j] = toolbox.rand_h1()
			elif key == 'h2':
				ind[key][j] = toolbox.rand_h2()
			elif key == 'linker_rots':
				ind[key][j] = toolbox.rand_angle()
			elif key == 'mutamers':
				ind[key][j] = toolbox.rand_aa()
			elif key == 'rotamers':
				ind[key][j] = toolbox.rand_rotamer()		
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
args = parser.parse_args()
wd = os.path.dirname(os.path.realpath(sys.argv[0]))

# Get Chimera params
residues = chimera.selection.savedSels['mutable'].residues()
base_at = chimera.selection.savedSels['base'].atoms()[0]
mol = base_at.molecule
ligand = base_at.residue
anchor = chimera.selection.savedSels['anchor'].atoms()[0]
base_geom = chimera.idatm.typeInfo[anchor.idatmType].geometry - 1
cbase = [base_at] + list(hyde5.atomsBetween(base_at, anchor)) + [anchor]

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
deap.creator.create("FitnessMax", deap.base.Fitness, weights=(1.0, -1.0, -1.0))
deap.creator.create("Individual", dict, fitness=deap.creator.FitnessMax)

# Operators
toolbox = deap.base.Toolbox()
toolbox.register("rand_angle", random.randint, 0, 359)
toolbox.register("rand_h1", random.randint, 0, base_geom-1)
toolbox.register("rand_h2", random.randint, 0, 2)
toolbox.register("rand_aa", random.randint, 0, len(aminoacids)-1)
toolbox.register("rand_rotamer", random.randint, 0, 8)
toolbox.register("rand_linker", random.randint, 0, len(linkers)-1)
toolbox.register("rand_fragment", random.randint, 0, len(fragments)-1)

# Genes
toolbox.register("molecule", deap.tools.initCycle, list,
	[toolbox.rand_linker, toolbox.rand_fragment], n=1)
toolbox.register("h1", deap.tools.initRepeat, list,
	toolbox.rand_h1, n=1)
toolbox.register("h2", deap.tools.initRepeat, list,
	toolbox.rand_h2, n=1)
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
	toolbox.molecule, toolbox.h1, toolbox.h2, toolbox.linker_rots,
	toolbox.mutamers, toolbox.rotamers)
toolbox.register("population", deap.tools.initRepeat, 
	list, toolbox.individual)

# Aliases for algorithm
toolbox.register("evaluate", evalCoord)
toolbox.register("mate", hetCxOnePoint)
toolbox.register("mutate", hetMutation, indpb=0.05)
toolbox.register("select", deap.tools.selNSGA2)

def main():
	pop = toolbox.population(n=args.pop)
	hof = deap.tools.HallOfFame(1)
	stats = deap.tools.Statistics(lambda ind: ind.fitness.values[1])
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
	evalCoord(hof[0], close=False)
	print("Best individual is: %s\nwith fitness: %s" % (hof[0], hof[0].fitness))

	# test =  {'mutamers': [0, 12], 'h2': [1], 'h1': [0], 'molecule': [0, 0], 
	# 'rotamers': [5, 1], 'linker_rots': [85, 342, 103, 240, 171, 215, 328, 159]}
	# print "Individual:\n{0}\nFitness:\n{1}".format(test, evalCoord(test, close=False))
