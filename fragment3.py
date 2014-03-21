"""
Select a terminal atom and run:
runscript fragment.py <SMILES string>
"""

import chimera
import sys

def insertMol(mol2, target=None, join=True, inplace=True, alpha=120.0):
	from chimera import UserError

	if not target:
		try: 
			target = chimera.selection.currentAtoms()[0]
		except KeyError:
			raise UserError("Please, select an H atom")

	# if not target.element.name == 'H':
	# 	raise UserError("Selected target must be an H atom")

	tmpl = chimera.openModels.open(mol2, type="Mol2", hidden=False)[0]

	# Place it nicely
	if inplace:
		from chimera import Xform as x
		from AddAttr import addAttributes

		addAttributes(mol2.replace('.mol2','.attr'), raiseTool=False)

		for a in tmpl.atoms:
			if a.element.name == 'H':
				tmpl.deleteAtom(a)
				continue
			if a.anchor in (3,5,7): anchor = a 
			if a.anchor in (1,5,6): axis_start = a 
			if a.anchor in (2,7,8): axis_end = a 

		# align target+anchor
		dv = target.coord() - anchor.coord()  # translation vector
		t = x.translation(dv) # translation xform matrix
		for a in tmpl.atoms: 
			a.setCoord(t.apply(a.coord()))

		#discard H atom and set actual target
		new_target = target.neighbors[0]
		target.molecule.deleteAtom(target)
		target = new_target

		# rotate params
		zero = chimera.Point(0.0, 0.0, 0.0)
		bond_axis = target.coord() - anchor.coord()
		mol_axis = axis_end.coord() - axis_start.coord()
		delta = chimera.angle(mol_axis,bond_axis) - alpha
		axis = chimera.cross(mol_axis, bond_axis)
		if axis.data() == (0.0,0.0,0.0):
			axis = chimera.cross(mol_axis, bond_axis + chimera.Vector(1,0,0))

		try: #actual rotation
			r = x.translation(anchor.coord() - zero) # move to origin
			r.multiply(x.rotation(axis, delta)) # rotate
			r.multiply(x.translation(zero -  anchor.coord())) # return to orig pos
			for a in tmpl.atoms:
				a.setCoord(r.apply(a.coord()))

		except ValueError: #this means rot vector is null  (already parallel)
			pass

		if join:
			from chimera.molEdit import addAtom, addBond
			from chimera.misc import getPseudoBondGroup 
			# convert PseudoBonds to regular bonds
			pseudobonds = getPseudoBondGroup("coordination complexes of %s (%s)" % 
				(tmpl.name, tmpl), associateWith=[tmpl]).pseudoBonds
			if pseudobonds:
				for pb in pseudobonds:
					addBond(pb.atoms[0],pb.atoms[1])

			# rename atoms accordingly
			oldRes = target.residue
			index = getHighestAtomIndices(oldRes)
			for a in tmpl.atoms:
				try:
					index[a.element.name] += 1
					a.name = a.element.name + str(index[a.element.name])
				except KeyError:
					index[a.element.name] = 1
			built_atoms = []
			sprouts = [ anchor ] # start to grow from seed
			while sprouts:
				sprout = sprouts.pop(0) # get first atom
				if sprout.name in oldRes.atomsMap:
					#print sprout.name  + " already in"
					target = oldRes.atomsMap[sprout.name][-1]
				else:
					#print sprout.name  + " not in"
					built = addAtom(sprout.name, sprout.element, oldRes,
							sprout.coord(), None, bondedTo=target)
					built.anchor = sprout.anchor
					target = built
					built_atoms.append(built)

				for a in sprout.neighbors: # get neighbors (maybe sprout.neighbors works too)
					if a.name not in oldRes.atomsMap:
						needBuild = True
					else:
						# atom is already present, but it can be part of a cycle
						# if we get to it it's because another atom is linking it
						needBuild = False
						built = oldRes.atomsMap[a.name][-1]
					if needBuild:
						built = addAtom(a.name, a.element, oldRes,
							a.coord(), None, bondedTo=target)
						built.anchor = a.anchor
						built_atoms.append(built)
						# if a has more than one neighbor:
						if len(a.neighbors) > 1:
							sprouts.append(a) # this new atom can be a new sprout

					if built not in target.bondsMap: #link!
						addBond(target, built)	
			chimera.openModels.close([tmpl])
			return built_atoms
	return tmpl.atoms

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


#############################
if __name__ == '__main__':
	insertMol(sys.argv[1])
