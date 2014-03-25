## LEGO script

# The idea is to build a ligand molecule based on
# different building blocks such as:
# ANCHOR--[series of linkers]--[series of fragments]

import chimera
from chimera import selection as sel
import fragment3 as frag
import os, sys, random

def getBase():
	return sel.savedSels['base'].atoms(), sel.savedSels['anchor'].atoms()[0]

def getMol2Files(path):
	# default path sould be wd/mol2/<linkers | fragment>
	return [ os.path.join(path,fn) for fn in next(os.walk(path))[2] if fn.endswith('.mol2') ]

def clearBase(base, anchor, get_coords=False, xform=False):
	not_base = [ a for a in base[0].residue.atoms if a not in base ]
	if get_coords:
		if xform:
			coords = [ a.xformCoord() for a in not_base ]
		else:
			coords = [ a.coord() for a in not_base ]
			
	[ base[0].molecule.deleteAtom(a) for a in not_base ]
	
	if get_coords:
		return coords


#####
base, anchor = getBase()
wd = os.path.dirname(os.path.realpath(sys.argv[0]))
linkers = getMol2Files(wd +'/mol2/linkers/')
fragments = getMol2Files(wd + '/mol2/fragments/')

l = random.randint(0, len(linkers)-1)
f = random.randint(0, len(fragments)-1)
clearBase(base,anchor)
linker = frag.insertMol(linkers[l], target=anchor, join=True, inplace=True)
linker_anchor = [ a for a in linker if a.anchor in (4,6,8) ]
fragment = frag.insertMol(fragments[f], target=linker_anchor[0], alpha=-120.)
