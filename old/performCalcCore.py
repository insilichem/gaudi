# Chimera script

import FindHBond as fh
import DetectClash as fc
from chimera import replyobj

def performCalcCore(hbonds=True, clashes=True):
	## Calculates H bonds and clashes for current position of current selection
	## Using native Chimera-Python interface

	results = []
	if (hbonds):
		

		replyobj.status("HBONDS OK")
	if (clashes):
		pass

	return results	

performCalcCore()

def customFindHBonds() :

	fh.hbonds = findHBonds(cacheDA=True)