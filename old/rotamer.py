## Rotamer switch
# Core implementation of swapaa

def optimizeRotamer(residue,num=5,clashes=None):
	#TODO save current rotamer/residue
	import Rotamers as r
	#initial conditions
	if clashes == None:
		clashes = countClashes(residue.atoms)[1]
		print "Found {0} contacts.".format(clashes)
	if clashes == 0:
		print "Rotamer is OK"
		return
	print "Optimizing"
	rotamers = r.getRotamers(residue)[1][:num]
	
	bestRot = None
	for rotamer in rotamers:
		print "Trying with rotamer chi " + str(rotamer.chis)
		#create ghost model to test clashes
		chimera.openModels.add([rotamer],sameAs=residue.molecule,hidden=True)
		#r.useRotamer(residue,[rotamer])
		allClashes  = countClashes(rotamer.atoms,intraRes=False,interSubmodel=True)[0]
		currClashes = 0
		resAtoms = set(residue.atoms)
		for ra in rotamer.atoms:
			if ra.name in ("CA", "N", "CB"):
				# any clashes of CA/N/CB are already clashes of
				# base residue (and may mistakenly be thought
				# to clash with "bonded" atoms in nearby
				# residues
				continue
			if ra not in allClashes:
				continue
			for ca, clash in allClashes[ra].items():
				if ca in resAtoms:
					continue
				if ca.molecule.id != residue.molecule.id:
					continue
				currClashes +=1
		print currClashes
		if currClashes < clashes:
			print "This rotamer is better! "
			clashes = currClashes
			bestRot = rotamer
		#remove ghost
		chimera.openModels.remove(rotamer)

	if bestRot != None:
		print bestRot.chis
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

def countClashes(atoms=[],intraRes=True,interSubmodel=False):
	#calculate clashes
	import DetectClash as dc
	if not atoms: atoms = chimera.selection.currentResidues()[0].atoms
	clashes = dc.detectClash(atoms,intraRes=intraRes,interSubmodel=interSubmodel)
	num_of_clashes = 0
	if clashes:
		for clashList in clashes.values():
			num_of_clashes += len(clashList)

	return clashes, num_of_clashes/2
#################

import chimera
from chimera import selection as sel

res = sel.currentResidues()[0]
best_rot = optimizeRotamer(res, 8)