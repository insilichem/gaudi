###########################################################
# CHIMERA SCRIPT
# Rotate selected bonds, find clashes and H bonds
# By Jaime RG <jaime.rogue@gmail.com>, in UAB
###########################################################

## VERSION 4a
## drop command line input, use core python
# [X] Rotation
# [ ] FindHbonds
# [ ] FindClashes


## KNOWN ISSUES

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

def getBonds(bonds=[], select=True, midas=False):
	## Parse current selection and write Midas-compatible selection string
	## using :id@first-atom@second-atom	
	results = []
	if midas: 
		if select: 
			for bond in selection.currentBonds():
				results.append(":" + str(bond.atoms[0].residue.id) + "@" + 
					str(bond.atoms[0].name) + "@" + str(bond.atoms[1].name))
		if bonds: 
			for bond in bonds:
				results.append(":" + str(bond.atoms[0].residue.id) + "@" + 
					str(bond.atoms[0].name) + "@" + str(bond.atoms[1].name))
		return results
	if select:
		return selection.currentBonds()
	else:
		raise "getBonds needs either a selection or a set of bonds"

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

def getReplyLog (chunk=100, flush=False):
	from chimera import dialogs, replyobj, UserError

	if "IDLE" in replyobj.currentReply().__doc__:
		raise UserError("Sorry, IDLE must be closed.")

	#dialogs.find('reply',create=1).enter()
	r = dialogs.find('reply').text.get('1.0','end')[-chunk:].split("\n")
	if flush:
		replyobj.currentReply().flush()
	return r


def performCalc():
	## Calculates H bonds and clashes for current position of current selection
	if chimera.nogui: #read from stdout
		# Capture STDOUT (somehow)
		from cStringIO import StringIO
		backup = sys.stdout
		sys.stdout = StringIO() #capture output
		#calculate hbonds & clashes
		rc("findhbond selRestrict any cacheDA true log true")
		hbonds_output = sys.stdout.getvalue()
		sys.stdout.close()
		sys.stdout = StringIO()
		rc("findclash ligand makePseudobonds false log true")
		clashes_output = sys.stdout.getvalue()
		sys.stdout.close()
		sys.stdout = backup # restore original stdout
		# Release STDOUT

		# Parse STDOUT
		clashes = clashes_output.splitlines()[7].split()[0]
		hbonds = 0
		for line in hbonds_output.splitlines():
			if line[:1] == ':' :
				hbonds += 1
	else: #read from reply log
		#calculate H bonds
		rc("findhbond selRestrict any cacheDA true log true")
		hbonds = getReplyLog(50,flush=True)[-3].split()[0]
		#calculate clashes
		rc("findclash ligand makePseudobonds false log true")
		clashes = getReplyLog(50,flush=True)[-3].split()[0]
	
	if clashes == "No" :
		clashes = 0
	
	# return hbonds, clashes
	return hbonds, clashes

def performCalcCore(hbonds=True, clashes=True):
	## Calculates H bonds and clashes for current position of current selection
	## Using native Chimera-Python interface
	results = []
	if (hbonds):
		pass
	if (clashes):
		pass

	return results	

def takeSnapshot(indices,wd):
	# takes a photo of current position
	rc("copy file " + wd + "img/" + "_".join(map(str,indices)) + ".jpg")

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
selected_bonds = getBonds()
anchor = selection.currentAtoms()[0]
for bond in selected_bonds:
	createRotation(bond,anchor)

from BondRotMgr import bondRotMgr
for rot in bondRotMgr.rotations.values():
	print rot, rot.anchorSide

# Select ligand that contains first bond
rc("sel :" + str(selected_bonds[0].atoms[0].residue.id))
rc("namesel ligand") #save sel
rc("focus ligand")

output = open(wd + "data.csv", 'w')
output.write("#Chosen bonds;--->;" + 
	" & ".join(getBonds(selected_bonds,midas=True)) + "\n")
output.write("Coordinates;Contacts;H bonds;\n")

## CORE
# Explore and evaluate

i = 0
for prev, curr in parComparison(product(xrange(0,360,args.step), 
										repeat=len(selected_bonds))):
	print "----\nCoordinate " + str(curr)
	# Calculate starting position ALWAYS
	if (prev == None) :
		clashes, hbonds = performCalc()
		# if (int(clashes) == 0): 
		# 	takeSnapshot(curr, wd)
		output.write(str(curr) + ";" + str(clashes) + ";" + str(hbonds) + "\n")
		continue

	for a in xrange(0,len(prev)):
		if (prev[a] != curr[a]):
			# rotate
		 	rotate(selected_bonds[a], curr[a] - prev[a])

 	# get score of current position
	clashes, hbonds = performCalc()
	# if (int(clashes) == 0): 
	#	takeSnapshot(curr, wd)
	output.write(str(curr) + ";" + str(clashes) + ";" + str(hbonds) + "\n")

	# progress
	i += 1
	total_rotations = (360/args.step)**len(selected_bonds)
	progress = 100*i/total_rotations
	if not progress % 5:
		print "Progress: " + str(progress) + "% done.\n"

output.close()

#rc("~sel")
replyobj.status("OK. Check results in data.csv in script dir", 
	blankAfter=5, log=True)