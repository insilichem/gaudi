###########################################################
# CHIMERA SCRIPT
# Rotate selected bonds, find clashes and H bonds
# By Jaime RG <jaime.rogue@gmail.com>, in UAB
###########################################################

## VERSION 4b
# drop command line input, use core python
# [X] Rotation
# [X] FindHbonds
# [X] FindClashes

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

def getMidasBonds(bonds=[]):
	## Parse current selection and write Midas-compatible selection string
	## using :id@first-atom@second-atom	
	results = [] 
	if bonds: 
		for bond in bonds:
			results.append(":{0}@{1}@{2}".format(bond.atoms[0].residue.id,
				bond.atoms[0].name, bond.atoms[1].name))
	else: # selection mode 
		for bond in selection.currentBonds():
			results.append(":{0}@{1}@{2}".format(bond.atoms[0].residue.id,
				bond.atoms[0].name, bond.atoms[1].name))
	return results

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

def clearRotation(bonds, allbonds=True):
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

	if anchor in br.atoms:
		br.anchorSide = anchor
	else: 
		br.anchorSide = findNearest(anchor, br.atoms)

def rotate(bond, degrees, absolute=False):
	## Python wrapper for rotation Chimera command
	from BondRotMgr import bondRotMgr
	br = bondRotMgr.rotationForBond(bond, create=False)

	if absolute:
		br.set(degrees)
	else:
		br.increment(degrees)

def countHBonds(model=chimera.openModels.list(modelTypes=[chimera.Molecule]),
		selRestrict=True):
	## Calculates H bonds and clashes for current position of current selection
	
	#calculate H bonds
	import FindHBond as fhb
	hbonds = fhb.findHBonds(model, cacheDA=True)
	if selRestrict:
		hbonds = fhb.base._filterBySel(
			hbonds,chimera.selection.currentAtoms(asDict=True), 'any')
	
	return len(hbonds)

def countClashes(atoms=chimera.selection.currentAtoms()):
	#calculate clashes
	import DetectClash as dc
	clashes = dc.detectClash(atoms,intraRes=True)
	num_of_clashes = 0
	if clashes:
		for clashList in clashes.values():
			num_of_clashes += len(clashList)

	return num_of_clashes/2


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
	ends = set([atom1, atom2])
	return chain|ends

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

def getRotamers(residue, num):
	import Rotamers as r
	rotamers = r.getRotamers(residue)[1]

	return rotamers[num:]

def useBestRotamer(residue,rotamers,clashes=0):
	import Rotamers as r
	#initial conditions
	#clashes = countClashes(residue.atoms)
	bestRotamer = residue

	#try rotamers
	for rotamer in rotamers:
		r.useRotamer(residue,[rotamer])
		if countClashes(residue.atoms) < clashes:
			clashes = countClashes(residue.atoms)
			bestRotamer = rotamer

	return bestRotamer


############################
#/ END OF FUNCTIONS
############################

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

# Activate rotation for selected bonds
selection.saveSel("initial")
selected_bonds = selection.currentBonds()
anchor = selection.currentAtoms()[0]
for bond in selected_bonds:
	createRotation(bond,anchor)

# Select current Residues
selection.setCurrent(selection.currentResidues())
selection.saveSel("ligand") #save sel

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
		hbonds = countHBonds()
		clashes = countClashes()
		coord = '-'.join(map(str,curr))
		output.write("{0};{1};{2};\n".format(coord,clashes,hbonds))
		continue

	for a in xrange(0,len(prev)):
		if (prev[a] != curr[a]):
		 	rotate(selected_bonds[a], curr[a] - prev[a])

 	# get score of current position
	hbonds = countHBonds()
	clashes = countClashes()
	coord = '-'.join(map(str,curr))

	output.write("{0};{1};{2}\n".format(coord,clashes,hbonds))

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