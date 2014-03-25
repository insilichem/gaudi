"""
Select a terminal atom and run:
runscript fragment.py <SMILES string>
"""

import chimera
import sys

def insertMol(mol2=None, pubchemid=None, smiles=None, 
	target=None, join=True, addlinker=True, inplace = True):
	import BuildStructure as bs	
	from chimera import match

	if not target:
		target = chimera.selection.currentAtoms()[0]
	model = target.molecule

	if mol2: #open mol2 file
		mol = chimera.openModels.open(mol2, type="Mol2", sameAs=model, hidden=True)[0]
	elif smiles:
		from BuildStructure import Smiles
		if join:
			smiles = "C"+smiles
		mol = Smiles.smiles2mol(smiles,resName=target.residue.type)
		chimera.openModels.add([mol],hidden=True,sameAs=model,)
	elif pubchemid:
		#from BuildStructure import PubChem
		pass
	for a in mol.atoms:
		if a.element.name == 'H':
			mol.deleteAtom(a)

	# Place it nicely
	if inplace:
		place(mol, target, 2., 0)

	return mol

def place(mol, target, dist, ang):
	# from http://plato.cgl.ucsf.edu/pipermail/chimera-dev/2010/000701.html

	from chimera import cross, angle, Point, Xform

	# fca = central fragment atom (it will be linked!)
	# fsa = fragment substituent atom (for angle positioning)
	# oa = atom from target (for distance)
	# mol = fragment residue
	# fca = mol.atoms[0]
	fca = list(terminii(mol))[0]
	fsa = target
	oa = target.neighbors[0]

	# make fca-oa distance = $dist
	dv = fca.coord() - oa.coord()
	dv.length = dist - dv.length
	for a in mol.atoms:
		a.setCoord(a.coord() + dv)

	# make fsa-fca-oa angle = $ang
	fscrd, fccrd, ocrd = [a.coord() for a in (fsa, fca, oa)]
	axis = cross(fscrd-fccrd, fccrd-ocrd) # axis the rotation will be around
	delta = ang - angle(fscrd, fccrd, ocrd)
	# rotation is about origin, so move center atom there
	v = fccrd - Point(0.0, 0.0, 0.0) 
	# translation to origin
	trans1 = Xform.translation(v) 
	v.negate() ## ????
	# translation back from origin
	trans2 = Xform.translation(v) 
	trans1.multiply(Xform.rotation(axis, delta))
	trans1.multiply(trans2)
	for a in mol.atoms:
		a.setCoord(trans1.apply(a.coord()))

def buildSprout(template, target, close=True):
	from chimera.molEdit import addAtom, addBond
	oldRes = target.residue
	# seed = list(terminii(mol))[0]
	seed = template.atoms[0]

	index = getHighestAtomIndices(oldRes)
	for a in template.atoms:
		index[a.element.name] += 1
		a.name = a.element.name + str(index[a.element.name])

	builtSprout = [ addAtom(seed.name, seed.element, 
			oldRes, seed.coord(), None, 
			bondedTo=target) ]
	
	# TODO: Instead of following the list of atoms
	# we could try to follow the whole tree, branch by branch
	for sprout in template.atoms:
		builtSprout.append(oldRes.atomsMap[sprout.name][-1])
		#builtSprout.append(new_sprout)
		for a, b in sprout.bondsMap.items():
			if a.name not in [ x.name for x in builtSprout[-1].neighbors ]: 
				built_a = addAtom(a.name, a.element, oldRes,
					a.coord(), None, bondedTo=builtSprout[-1])
				
				if built_a not in builtSprout[-1].bondsMap:
					addBond(builtSprout, built_a)
	if close:
		chimera.openModels.close([template])

def getHighestAtomIndices(r):
	results = {}
	for a in r.atoms:
		if a.name[1:].isdigit():
			atom = a.name[:1]
			num = int(a.name[1:])
			if atom not in results:
				results[atom] = num
				
			elif results[atom] < num:
				results[atom] = num
				
		elif a.name[2:].isdigit():
			atom = a.name[:2]
			num = int(a.name[2:])
			if atom not in results:
				results[atom] = num
				
			elif results[atom] < num:
				results[atom] = num
	return results

def terminii(residue):
	# returns single-bond atoms
	for at in residue.atoms:
		if at.numBonds == 1:
			yield at

target = chimera.selection.currentAtoms()[0]

mol = insertMol(mol2=sys.argv[1], target=target, inplace=True)
buildSprout(mol, target)º