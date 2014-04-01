###########################################################
# CHIMERA SCRIPT
# Rotate selected bonds, find clashes and H bonds
# By Jaime RG <jaime.rogue@gmail.com>, in UAB
###########################################################

## TODO
# 1. Save and Restore original state
# 2. Rebuild given coodinates from scratch
# 3. BUILD COORDINATES MAP
# 	 Maybe we don't want a fixed step?
#######

## KNOWN ISSUES
# 1. File output and parsing is slow
#
#######

import os
import sys
from itertools import tee, islice, chain, izip, product
import argparse
from chimera import runCommand as rc #cmd line features
from chimera import replyobj

def getBonds():
	## Parse current selection and write Midas-compatible selection string
	## using :id@first-atom@second-atom 
	results = []
	for bond in chimera.selection.currentBonds():
			results.append(":" + str(bond.atoms[0].residue.id) + "@" + str(bond.atoms[0].name) + "@" + str(bond.atoms[1].name))
	return results

def parComparison(some_iterable):
	## Returns the iterable in tuples of two, lazily
    prevs, items = tee(some_iterable, 2)
    prevs = chain([None], prevs)
    return izip(prevs, items)

def activateRotation(selected_bonds, reverse=False):
	## Activates selected bonds for rotation
	## TODO: Smart reverse option, based on picked anchor
	i = 0
	for bond in selected_bonds:
		try:
			if (reverse):
				rc("rotation " + str(i) + " reverse " + bond)
			else:
				rc("rotation " + str(i) + " " + bond)
		except:# chimera.UserError:
			rc("~rotation " + str(i))
			if (reverse):
				rc("rotation " + str(i) + " reverse " + bond)
			else:
				rc("rotation " + str(i) + " " + bond)
		i += 1

def rotate(index, degrees):
	## Python wrapper for rotation Chimera command
	rc("rotation " + str(index) + " " + str(degrees))

def performCalc(indices, wd):
	## Calculates H bonds and clashes for current position of current selection
	
	#calculate H bonds
	fn_hbonds = wd + "tmp/degree__" + str(indices) + ".hbonds.txt"
	rc("findhbond selRestrict any cacheDA true saveFile \"" + fn_hbonds + "\"")
	#calculate clashes
	fn_clash = wd + "tmp/degree__" + str(indices) + ".clash.txt"
	rc("findclash ligand makePseudobonds false saveFile \"" + fn_clash + "\"")

	#parse clashes
	f1 = open(fn_clash, 'r')
	clashes = f1.readlines()
	num_of_clashes = clashes[6].split(" ", 1)
	f1.close()
	#parse H bonds; careful! the source needs to be modified for this to work
	#since it will look for a summary line not present in original sc
	f2 = open(fn_hbonds, 'r')
	hbonds = f2.readlines()
	num_of_hbonds = hbonds[0].split(" ", 1)
	f2.close()

	return num_of_clashes[0], num_of_hbonds[0]

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

###########################
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


#working directory
wd = os.path.dirname(os.path.realpath(sys.argv[0])).replace('\\', '/') + '/'
rc("cd " + wd)

# Check bond is selected
try:
	chimera.selection.currentBonds()[0]
except IndexError:
	raise IndexError("Select one or more bonds.")
	replyobj.error("Select one or more bonds.", color="red")

# Activate rotation for selected bonds
selected_bonds = getBonds()
activateRotation(selected_bonds, args.reverse)


# Select ligand that contains bond
rc("sel :" + str(chimera.selection.currentBonds()[0].atoms[0].residue.id))
rc("namesel ligand") #save sel
# Focus
rc("focus ligand")

# Prepare wd and files
if not os.path.isdir(wd + "tmp/"):
	os.makedirs(wd + "tmp/")
if not os.path.isdir(wd + "img/"):
	os.makedirs(wd + "img/")
else: # remove previous images if any
	[ os.remove(wd + "img/" + f) for f in os.listdir(wd + "img/") ]

# Results file
output = open(wd + "data.csv", 'w')
output.write("Coordinates;Contacts;H bonds;Chosen bonds ->;" + ";".join(selected_bonds) + "\n")

## CORE
# Explore and evaluate
i = 0

for prev, curr in parComparison(product(xrange(0,360,args.step), repeat=len(selected_bonds))):
	print "----\nCoordinate " + str(curr)
	# Calculate starting position ALWAYS
	if (prev == None) :
		clashes, hbonds = performCalc(curr,wd)
		# if (int(clashes) == 0): 
		# 	takeSnapshot(curr, wd)
		output.write(str(curr) + ";" + str(clashes) + ";" + str(hbonds) + "\n")
		continue

	for a in xrange(0,len(prev)):
		if (prev[a] != curr[a]):
			# rotate
		 	rotate(a, curr[a] - prev[a])

 	# get score of current position
	clashes, hbonds = performCalc(curr,wd)
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

# clean tmp and clear selection
[ os.remove(wd + "tmp/" + f) for f in os.listdir(wd + "tmp/") ]
os.rmdir(wd + "tmp/")
rc("~sel")


replyobj.status("OK. Check results in data.csv in script dir", blankAfter=5, log=True)