"""
Select a terminal atom and run:
runscript fragment.py <SMILES string>
"""

import chimera
import sys

def insertMol(mol2=None, pubchemid=None, smiles=None, 
	target=None, join=True, inplace = False, attr=True):
	from chimera import UserError
	from AddAttr import addAttributes
	if not target:
		target = chimera.selection.currentAtoms()[0]
	if not target.element.name == 'H':
		raise UserError("Selected target must be an H atom")

	model = target.molecule

	if mol2: #open mol2 file
		mol = chimera.openModels.open(mol2, type="Mol2", hidden=False)[0]
		if attr: #contains axis data
			addAttributes(mol2.replace('mol2','attr'), raiseTool=False)
	elif smiles:
		from BuildStructure import Smiles
		if join:
			smiles = "C"+smiles
		mol = Smiles.smiles2mol(smiles,resName=target.residue.type)
		chimera.openModels.add([mol],hidden=True)
	elif pubchemid:
		#from BuildStructure import PubChem
		pass
	
	for a in mol.atoms:
		if a.element.name == 'H':
			mol.deleteAtom(a)
		# custom attributes!
		elif a.anchor == 1:
			anchor = a
		elif a.anchor == 2:
			anchor_end = a

	# Place it nicely
	if inplace:
		from chimera import Xform as x

		#discard H atom, but save relevant info
		new_target = target.neighbors[0]
		insertion_coords = target.coord()
		target.molecule.deleteAtom(target)
		target = new_target

		# align target+anchor
		dv = insertion_coords - anchor.coord()  # translation vector
		t = x.translation(dv) # translation xform matrix
		for a in mol.atoms:
			a.setCoord(t.apply(a.coord()))

		# rotate params
		alpha = 180.0 # it should depend on target atom geometry
		zero = chimera.Point(0.0, 0.0, 0.0)
		bond_axis = target.coord() - anchor.coord()
		mol_axis = anchor_end.coord() - anchor.coord()
		delta = alpha - chimera.angle(mol_axis,bond_axis)
		axis = chimera.cross(mol_axis, bond_axis)

		try: #actual rotation
			r = x.translation(anchor.coord() - zero) # move to origin
			r.multiply(x.rotation(axis, - delta)) # rotate
			r.multiply(x.translation(zero -  anchor.coord())) # return to orig pos
		
			for a in mol.atoms:
				a.setCoord(r.apply(a.coord()))
		except ValueError: #this means cross vector is null  (already parallel)
			pass
	return mol, target

def buildSprout(template, target, close=True):
	from chimera.molEdit import addAtom, addBond

	# rename atoms accordingly
	oldRes = target.residue
	index = getHighestAtomIndices(oldRes)
	for a in template.atoms:
		index[a.element.name] += 1
		a.name = a.element.name + str(index[a.element.name])


	# try to get a terminal atom that acts as seed for growth
	try :
		seed = list(terminii(template))[0]
	except IndexError:
		seed = template.atoms[0]

	sprouts = [ seed ] # start to grow from seed
	while sprouts:
		sprout = sprouts.pop(0) # get first atom
		if sprout.name in oldRes.atomsMap:
			print sprout.name  + " already in"
			target = oldRes.atomsMap[sprout.name][-1]
		else:
			print sprout.name  + " not in"
			built = addAtom(sprout.name, sprout.element, oldRes,
					sprout.coord(), None, bondedTo=target)
			target = built

		for a in sprout.neighbors: # get neighbors (maybe sprout.neighbors works too)
			if a.name not in oldRes.atomsMap:
				needBuild = True
			else:
				# atom is already present, but it can be part of a cycle
				# if we get to it it's because another atom is linking it
				needBuild = False
				built_a = oldRes.atomsMap[a.name][-1]
			if needBuild:
				built_a = addAtom(a.name, a.element, oldRes,
					a.coord(), None, bondedTo=target)
				# if a has more than one neighbor:
				if len(a.neighbors) > 1:
					sprouts.append(a) # this new atom can be a new sprout

			if built_a not in target.bondsMap: #link!
				addBond(target, built_a)	

	if close:
		chimera.openModels.close([template])

def getHighestAtomIndices(r):
	results = { 'C': 0, 'H': 0, 'N' : 0, 'O': 0 }
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


def main(mol2, target):
	mol, new_target = insertMol(mol2=mol2, target=target, inplace=True)

	buildSprout(mol, new_target)

target = chimera.selection.currentAtoms()[0]
main(sys.argv[1], target)
