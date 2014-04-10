"""
Select a terminal atom and run:
runscript fragment.py <SMILES string>
"""

import chimera
import sys

def insertMol(mol2, target=None, join=True, inplace=True,
	p2b=True, alpha=None, alpha2=None, alpha3=None, dihedral=None):
	from chimera import UserError

	if not target:
		try: 
			target = chimera.selection.currentAtoms()[0]
		except KeyError:
			raise UserError("Please, select an H atom")

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

		#discard H atom and set actual target
		if target.element.number == 1:
			H_coord = target.coord()
			new_target = target.neighbors[0]
			target.molecule.deleteAtom(target)
			target = new_target
		else:
			from BuildStructure import changeAtom
			try:
				geometry = target.getIdatmInfoMap()[target.idatmType].geometry
			except KeyError:
				geometry = 4
			try:
				target, newH = changeAtom(target, target.element, geometry, 
					target.numBonds + 1) # TODO should be `geometry - target.numBonds`
				# target = newAtom[0]
				# newHs = newAtom[1:] # TODO consider all possible H -> newAtom[1:]
				# newH = newAtom[1]
			except IndexError:
				newH = [a for a in target.neighbors if a.numBonds == 1][-1]
				print "Warning, arbitrary H"
			
			# fix bond length
			from BuildStructure import elementRadius
			dv = newH.coord() - target.coord()
			dv.length = elementRadius[anchor.element] + elementRadius[target.element]
			diff = dv - (newH.coord() - target.coord())
			newH.setCoord(newH.coord() + diff)
			H_coord = newH.coord()
			target.molecule.deleteAtom(newH)
			
		# align target+anchor
		dv = H_coord - anchor.coord()  # translation vector
		t = x.translation(dv) # translation xform matrix
		for a in tmpl.atoms: 
			a.setCoord(t.apply(a.coord()))

		# # rotate params for angle
		zero = chimera.Point(0.0, 0.0, 0.0)

		# Critical atoms for rotation
		d5 = anchor.neighbors[0]
		d4 = anchor
		d3 = target
		d2 = [ a for a in target.neighbors if a.element.number != 1 ][0]
		d1 = [ a for a in d2.neighbors if a != target ][0]

		if alpha:
			# rotate params for angle
			axis_a = d2.coord() - d3.coord()
			axis_b = d4.coord() - d3.coord()
			delta = chimera.angle(d2.coord(),d3.coord(),d4.coord()) - alpha
			axis = chimera.cross(axis_a, axis_b)
			if axis.data() == (0.0,0.0,0.0):
				axis = chimera.cross(axis_a, axis_b + chimera.Vector(1,0,0))
				print "Warning, had to choose arbitrary normal vector"
			try: #actual rotation
				r = x.translation(target.coord() - zero) # move to origin
				r.multiply(x.rotation(axis, - delta)) # rotate
				r.multiply(x.translation(zero -  target.coord())) # return to orig pos
				for a in tmpl.atoms:
					a.setCoord(r.apply(a.coord()))
			except ValueError: #this means rot vector is null  (already parallel)
				raise ValueError("Rotation 1 failed. Normal vector is null")

		if alpha3:
			# rotate params for angle
			axis_a = d3.coord() - d4.coord()
			axis_b = axis_end.coord() - d4.coord()
			delta = chimera.angle(d3.coord(),d4.coord(),axis_end.coord()) - alpha3
			axis = chimera.cross(axis_a, axis_b)
			if axis.data() == (0.0,0.0,0.0):
				axis = chimera.cross(axis_a, axis_b + chimera.Vector(1,0,0))
				print "Warning, had to choose arbitrary normal vector"
			try: #actual rotation
				r = x.translation(d4.coord() - zero) # move to origin
				r.multiply(x.rotation(axis, - delta)) # rotate
				r.multiply(x.translation(zero -  d4.coord())) # return to orig pos
				for a in tmpl.atoms:
					a.setCoord(r.apply(a.coord()))
			except ValueError: #this means rot vector is null  (already parallel)
				raise ValueError("Rotation 3 failed. Normal vector is null")
		if alpha2:
			axis_a = d3.coord() - d4.coord()
			axis_b = d5.coord() - d4.coord()
			delta = chimera.angle(d3.coord(),d4.coord(),d5.coord()) - alpha2
			axis = chimera.cross(axis_a, axis_b)
			if axis.data() == (0.0,0.0,0.0):
				axis = chimera.cross(axis_a, axis_b + chimera.Vector(1,0,0))
				print "Warning, had to choose arbitrary normal vector"
			try: #actual rotation
				r = x.translation(anchor.coord() - zero) # move to origin
				r.multiply(x.rotation(axis, - delta)) # rotate
				r.multiply(x.translation(zero - anchor.coord())) # return to orig pos
				for a in tmpl.atoms:
					a.setCoord(r.apply(a.coord()))
			except ValueError: #this means rot vector is null  (already parallel)
				raise ValueError("Rotation 2 failed. Normal vector is null")

		if dihedral: # rotate params for dihedral
			axis = d3.coord() - d2.coord()
			delta = chimera.dihedral(d1.coord(), d2.coord(), d3.coord(), d4.coord()) - dihedral
			try: #actual rotation
				r = x.translation(target.coord() - zero) # move to origin
				r.multiply(x.rotation(axis, - delta)) # rotate
				r.multiply(x.translation(zero - target.coord())) # return to orig pos
				for a in tmpl.atoms:
					a.setCoord(r.apply(a.coord()))
			except ValueError: #this means rot vector is null  (already parallel)
				raise ValueError("Dihedral rotation failed. Normal vector is null")

		if join:
			from chimera.molEdit import addAtom, addBond
			if p2b: # convert PseudoBonds to regular bonds
				pbgroup = chimera.misc.getPseudoBondGroup(
					"coordination complexes of %s (%s)" % 
					(tmpl.name, tmpl), associateWith=[tmpl])
				if pbgroup.pseudoBonds:
					for pb in pbgroup.pseudoBonds:
						addBond(pb.atoms[0],pb.atoms[1])
					pbm = tmpl.pseudoBondMgr()
					pbm.deletePseudoBondGroup(pbgroup)

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
					if a.element.number == 1: continue
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
	insertMol(sys.argv[1], target=chimera.selection.currentAtoms()[0], alpha2=120,
		alpha3= 180.)
