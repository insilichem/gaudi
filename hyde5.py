###########################################################
# CHIMERA SCRIPT
# Rotate selected bonds, find clashes and H bonds
# By Jaime RG <jaime.rogue@gmail.com>, in UAB
###########################################################

## VERSION 5a
# include rotamers

## KNOWN ISSUES
# !! Reported HBONDS in 4b < HBONDS in 4a /!\

## TODO
""" 1. Save and Restore original state
	2. Rebuild given coodinates from scratch
 	3. BUILD COORDINATES MAP
 		Maybe we don't want a fixed step?
"""

import os
import sys
from itertools import tee, islice, chain, izip, product
import argparse
import chimera
from chimera import runCommand as rc #cmd line features
from chimera import replyobj, UserError, selection

def getMidasBonds(bonds=None):
	## Parse current selection and write Midas-compatible selection string
	## using :id@first-atom@second-atom	
	# results = [] 
	if bonds: 
		for bond in bonds:
			yield ":{0}@{1}@{2}".format(bond.atoms[0].residue.id,
				bond.atoms[0].name, bond.atoms[1].name)
	else: # selection mode 
		for bond in selection.currentBonds():
			yield ":{0}@{1}@{2}".format(bond.atoms[0].residue.id,
				bond.atoms[0].name, bond.atoms[1].name)
	

def parComparison(some_iterable):
	## Returns the iterable in tuples of two, lazily
    prevs, items = tee(some_iterable, 2)
    prevs = chain([None], prevs)
    return izip(prevs, items)

def dialogMgr(name, close=False):
	# Close requested dialog if open
	from chimera import dialogs
	
	# for d in dialogs.activeDialogs():
	# 	if d.name == name:
	# 		d.Close()
	# 		break
	if close and dialogs.find(name): dialogs.find(name).Close()
	else: 
		dialogs.find(name, create = 1).enter()

def clearRotation(bonds=None, allbonds=True):
	## Activates selected bonds for rotation
	# bonds: list of chimera bond objects
	from BondRotMgr import bondRotMgr

	if allbonds and bondRotMgr.rotations:
		for b in bondRotMgr.rotations.values():
			b.destroy()
	elif bonds and bondRotMgr.rotations:
		for b in bonds:
			b.destroy()

def createRotation(bond, anchor):
	from BondRotMgr import bondRotMgr
	br = bondRotMgr.rotationForBond(bond)

	if anchor not in br.atoms:
		anchor = findNearest(anchor, br.atoms)
	if br._BondRotation__anchorSide != anchor:
		br._BondRotation__anchorSide = anchor
		br.anchorSide = anchor
		bondRotMgr.triggers.activateTrigger(bondRotMgr.REVERSED, br)

def rotate(bond, degrees, absolute=False):
	## Python wrapper for rotation Chimera command
	from BondRotMgr import bondRotMgr
	br = bondRotMgr.rotationForBond(bond, create=False)

	if absolute:
		br.set(degrees)
	elif degrees:
		br.increment(degrees)

def bondrot(bond, anchor, delta):
	# Simpler than the two above. Avoids triggers though
	# Will have to ask about its implications
	from chimera import BondRot
	br = BondRot(bond)
	anchor = findNearest(anchor, bond.atoms)
	br.angle = (delta, anchor)
	br.destroy()
	
def countHBonds(model=None,	sel=None, selRestrict=True,cache=False):
	## Calculates H bonds and clashes for current position of current selection
	if not model:
		model = chimera.openModels.list(modelTypes=[chimera.Molecule])
	if not sel:
		sel = chimera.selection.currentAtoms()
	#calculate H bonds
	import FindHBond as fhb
	hbonds = fhb.findHBonds(model, cacheDA=cache)
	if selRestrict:
		hbonds = fhb.base._filterBySel(hbonds, sel, 'any')
	
	return hbonds

def countClashes(atoms=None, test='others', intraRes=True,interSubmodel=False):
	#calculate clashes
	import DetectClash as dc
	if not atoms: 
		atoms = chimera.selection.currentResidues()[0].atoms
	clashes = dc.detectClash(atoms, test=test, intraRes=intraRes,
		interSubmodel=interSubmodel)
	num_of_clashes = 0
	if clashes:
		for clashList in clashes.values():
			num_of_clashes += len(clashList)

	return clashes, num_of_clashes/2


def atomsBetween(atom1, atom2):
	# find all connected atoms between, given two connected ends
	chain1 = atom1.neighbors
	chain2 = atom2.neighbors
	i = 0
	while i < len(chain1):
		a1 = chain1[i]
		if atom2 not in a1.neighbors:
			chain1.extend( [ a for a in a1.neighbors if a not in chain1 ])
		i += 1
	j = 0
	while j < len(chain2):
		a2 = chain2[j]
		if atom1 not in a2.neighbors:
			chain2.extend( [ a for a in a2.neighbors if a not in chain2 ])
		j += 1
	
	chain = set(chain1) & set(chain2)
	#ends = set([atom1, atom2])
	return chain

def findNearest(anchor, atoms):
	# return atom with mininum through-bond distance
	nearest = atoms[0]
	minimum = 10000

	for a in atoms:
		distance = len(atomsBetween(anchor, a))
		if distance < minimum:
			minimum = distance
			nearest = a

	return nearest

def optimizeRotamer(residue,newRes=None,num=5,clashes=None,replace=False):
	#TODO save current rotamer/residue
	import Rotamers as r
	#initial conditions
	if clashes == None:
		clashes = countClashes(residue.atoms)[1]
	if clashes == 0:
		return


	# else
	rotamers = r.getRotamers(residue)[1][:num]
	
	bestRot = None
	for rotamer in rotamers:
		#create ghost model to test clashes
		chimera.openModels.add([rotamer],sameAs=residue.molecule,hidden=True)
		#r.useRotamer(residue,[rotamer])
		allClashes  = countClashes(rotamer.atoms,intraRes=False,interSubmodel=True)[0]
		currClashes = 0
		resAtoms = set(residue.atoms)
		#filter clashes
		for ra in rotamer.atoms:
			if ra.name in ("CA", "N", "CB"):
				continue
			if ra not in allClashes:
				continue
			for ca, clash in allClashes[ra].items():
				if ca in resAtoms:
					continue
				if ca.molecule.id != residue.molecule.id:
					continue
				currClashes +=1
		if currClashes < clashes:
			clashes = currClashes
			bestRot = rotamer
		#remove ghost
		chimera.openModels.remove(rotamer)

	if bestRot != None:
		print "Using new rotamer with chi: " + str(bestRot.chis)
		if replace:
			r.useRotamer(residue,[bestRot])
		return bestRot
	return

def allKeys(d):
	for k, v in d.iteritems():
		if isinstance(v, dict):
			for k_ in allKeys(v): 
				yield k_
			yield k
		else:
			yield k

################
def main():
# ARGUMENT PARSING
	parser = argparse.ArgumentParser()
	parser.add_argument('-d', '--degrees',
						required=False,
						type=int,
						default=90, 
						dest="step",
						metavar="<number of degrees>",
						help="Degrees for each rotation" )
	args = parser.parse_args()
	#### /ARGUMENT PARSING


	#working directory
	wd = os.path.dirname(os.path.realpath(sys.argv[0])).replace('\\', '/') + '/'
	rc("cd " + wd)

	# Check bond is selected
	if not selection.currentBonds():
		raise UserError("Please, select one or more bonds.")
	if not selection.currentAtoms() or len(selection.currentAtoms()) > 1:
		raise UserError("Please, select only one atom to act as anchor.")
	# if args.mutable:
	# 	if not 'mutable' in selection.savedSels or not selection.savedSels['mutable'].residues():
	# 		raise UserError("Please, select which residues should be explored " +
	# 			"for rotamers. Name that selection as 'mutable'.")

	# Activate rotation for selected bonds
	selection.saveSel("initial")
	selected_bonds = selection.currentBonds()
	anchor = selection.currentAtoms()[0]
	for bond in selected_bonds:
		createRotation(bond,anchor)

	# Select current Residues
	selection.setCurrent(selection.currentResidues())
	selection.saveSel("ligand") #save sel
	# mutable_residues = selection.savedSels['mutable'].residues()

	output = open(wd + "data.csv", 'w')
	output.write("#Chosen bonds;--->;{0}\n".format(" & ".join(
		getMidasBonds(selected_bonds))))
	output.write("Coordinates;Contacts;H bonds\n")

	# Explore and evaluate
	i = 0
	for prev, curr in parComparison(product(xrange(0,360,args.step), 
											repeat=len(selected_bonds))):
		print "----\nCoordinate " + str(curr)
		# Calculate starting position ALWAYS
		if (prev == None) :
			hbonds = len(countHBonds())
			clashes, num_of_clashes = countClashes()
			coord = '-'.join(map(str,curr))
			output.write("{0};{1};{2};\n".format(coord,num_of_clashes,hbonds))
			continue

		for a in xrange(0,len(prev)):
			if (prev[a] != curr[a]):
			 	rotate(selected_bonds[a], curr[a] - prev[a])

	 	# get score of current position
		hbonds = len(countHBonds())
		clashes, num_of_clashes = countClashes()
		coord = '-'.join(map(str,curr))
		print '{0} H Bonds and {1} clashes'.format(hbonds,num_of_clashes)
		# if num_of_clashes:
		# 	for mutable in mutable_residues:
		# 		if set(mutable.atoms).intersection(set(allKeys(clashes))):
		# 			newRotamer = optimizeRotamer(mutable, 5, clashes=num_of_clashes, replace=False)
		# 			if newRotamer:
		# 				print 'Changing ' + str(mutable)
		# 				coord += ',[{0}]'.format(newRotamer.name)


		output.write("{0};{1};{2}\n".format(coord,num_of_clashes,hbonds))

		# progress
		i += 1
		total_rotations = (360/args.step)**len(selected_bonds)
		progress = 100*i/total_rotations
		if not progress % 5:
			print "Progress: {0}% done.".format(progress)

	output.close()

	selection.setCurrent(selection.savedSels["initial"])
	replyobj.status("OK. Check results in data.csv in script dir", 
		blankAfter=5, log=True)

############################
#/ END OF FUNCTIONS
############################
if __name__ == '__main__':
	main()