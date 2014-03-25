###########################################################
# CHIMERA SCRIPT
# Rotate selected bonds, find clashes and H bonds
# By Jaime RG <jaime.rogue@gmail.com>, in UAB
###########################################################

## VERSION 4
# drop command line input, use core python

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
from chimera import replyobj, UserError

def getBonds(bonds=[], selection=True,midas=False):
	## Parse current selection and write Midas-compatible selection string
	## using :id@first-atom@second-atom	
	results = []
	if midas: 
		if selection: 
			for bond in chimera.selection.currentBonds():
				results.append(":" + str(bond.atoms[0].residue.id) + "@" + str(bond.atoms[0].name) + "@" + str(bond.atoms[1].name))
		if bonds: 
			for bond in bonds:
				results.append(":" + str(bond.atoms[0].residue.id) + "@" + str(bond.atoms[0].name) + "@" + str(bond.atoms[1].name))
		return results
	if selection:
		return chimera.selection.currentBonds()
	else:
		raise "getBonds needs either a selection or a set of bonds"

def parComparison(some_iterable):
	## Returns the iterable in tuples of two, lazily
    prevs, items = tee(some_iterable, 2)
    prevs = chain([None], prevs)
    return izip(prevs, items)

def closeDialog(name):
	# Close requested dialog if open
	from chimera import dialogs
	for d in dialogs.activeDialogs():
		if d.name == name:
			d.Close()
			break
	#dialogs.find('name').Close()

def clearRotation(bonds=[], allbonds=True):
	## Activates selected bonds for rotation
	# bonds: list of chimera bond objects
	## TODO: Smart reverse option, based on picked anchor
	from BondRotMgr import bondRotMgr

	if allbonds and bondRotMgr.rotations:
		for b in bondRotMgr.rotations.values():
			b.destroy()
	elif bonds and bondRotMgr.rotations:
		for b in bonds:
			b.destroy()


def rotate(bond, degrees, absolute=False):
	## Python wrapper for rotation Chimera command
	from BondRotMgr import bondRotMgr
	b = bondRotMgr.rotationForBond(bond)
	if absolute:
		b.set(degrees)
	else:
		b.increment(degrees)

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

###########################activateRotation(selected_bonds)

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
parser.add_argument('-r', '--reverse',
					required=False, 
					default=False,
					dest="reverse",
					action='store_true',
					help="Wether to move heavy or light half" )
args = parser.parse_args()
#### /ARGUMENT PARSING


# Check bond is selected
if not chimera.selection.currentBonds():
	raise UserError("Please, select one or more bonds.")


# Activate rotation for selected bonds
selected_bonds = getBonds()
#activateRotation(selected_bonds)
#closeDialog('build structure')
clearRotation()
# Select ligand that contains bond
rc("sel :" + str(chimera.selection.currentBonds()[0].atoms[0].residue.id))
print "#Chosen bonds;--->;" + " & ".join(getBonds(selected_bonds,midas=True))

from BondRotMgr import bondRotMgr
## CORE
# Explore and evaluate
for prev, curr in parComparison(product(xrange(0,360,args.step), repeat=len(selected_bonds))):
	print "----\nCoordinate " + str(curr)
	if prev == None:
		continue
	for a in xrange(0,len(prev)):
		if (prev[a] != curr[a]):
			# rotate
		 	rotate(selected_bonds[a], curr[a] - prev[a])
		 	#print bondRotMgr.rotationForBond(selected_bonds[a]).get()
