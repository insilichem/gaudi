#!/usr/bin/python

# MOF3D
# Multi-Objective Force-Field-Free Docking
# A docking module for UCSF Chimera
# Jaime RGP <https://bitbucket.org/jrgp> @ UAB, 2014

import chimera, BuildStructure, AddAttr, SplitMolecule
from chimera.molEdit import addAtom, addBond
from itertools import product
#custom
import move, utils

def insert(mol2, target=None, join=True, p2b=True, inplace=True,
	alpha=None, alpha2=None, alpha3=None, dihedral=None):
	
	''' Expected parameters.
	target=None, join=True, p2b=True, inplace=True,
	alpha=None, alpha2=None, alpha3=None, dihedral=None
	'''

	tmpl = chimera.openModels.open(mol2, type="Mol2", hidden=False)[0]
	AddAttr.addAttributes(mol2.replace('.mol2','.attr'), models=[tmpl],raiseTool=False)
	
	for p in ('alpha', 'alpha2', 'alpha3', 'dihedral'):
		if not locals()[p] and hasattr(tmpl, p):
			locals()[p] = getattr(tmpl, p)
	
	# Place it nicely
	if inplace:
		if not target:
			try: 
				target = chimera.selection.currentAtoms()[0]
			except KeyError:
				raise chimera.UserError("Please, select an H atom")

		# Read and parse anchor attributes
		for a in tmpl.atoms:
			if a.element.name == 'H':
				tmpl.deleteAtom(a)
				continue
			if a.anchor in (3,5,7): anchor = a 
			if a.anchor in (1,5,6): axis_start = a 
			if a.anchor in (2,7,8): axis_end = a 

		#discard H atom and set actual target
		if target.element.number == 1:
			targetH, target = target, target.neighbors[0]
		else: #add new H based on target atom geometry
			try:
				geometry = chimera.idatm.typeInfo[target.idatmType].geometry
			except KeyError:
				geometry = 4
			try:
				target, targetH = BuildStructure.changeAtom(
								target, target.element, geometry, 
								target.numBonds + 1)
			except IndexError: #unpacking error; pick any terminal atom available
				targetH = [a for a in target.neighbors if a.numBonds == 1][-1]
				print "Warning, arbitrary H"
			
		# fix bond length
		dv = targetH.coord() - target.coord()
		dv.length = chimera.Element.bondLength(anchor.element, target.element)
		diff = dv - (targetH.coord() - target.coord())
		targetH.setCoord(targetH.coord() + diff)
		H_coord = targetH.coord()
		target.molecule.deleteAtom(targetH)
			
		# align target+anchor
		move.translate(tmpl, anchor, H_coord)
		# Critical atoms for rotation
		d5 = anchor.neighbors[0]
		d4 = anchor
		d3 = target
		d2 = [ a for a in target.neighbors if a.element.number != 1 ][0]
		d1 = [ a for a in d2.neighbors if a != target ][0]

		if alpha:
			move.rotate(tmpl, [d2,d3,d4], alpha)
		if alpha3:
			move.rotate(tmpl, [d3, axis_start, axis_end], alpha3)
		if alpha2:
			move.rotate(tmpl, [d3,d4,d5], alpha2)
		if dihedral:
			move.rotate(tmpl, [d1, d2, d3, d4], dihedral)
		if join:
			if p2b: # convert PseudoBonds to regular bonds
				utils.box.pseudobond_to_bond(tmpl)
			# rename atoms accordingly
			old_res = target.residue
			index = utils.box.highest_atom_indices(old_res)
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
				if sprout.name in old_res.atomsMap:
					target = old_res.atomsMap[sprout.name][-1]
				else:
					built = addAtom(sprout.name, sprout.element, old_res,
							sprout.coord(), None, bondedTo=target)
					built.anchor = sprout.anchor
					target = built
					built_atoms.append(built)

				for a in sprout.neighbors:
					if a.element.number == 1: continue
					if a.name not in old_res.atomsMap:
						needBuild = True
					else:
						# atom is already present, but it can be part of a cycle
						# if we get to it it's because another atom is linking it
						needBuild = False
						built = old_res.atomsMap[a.name][-1]
					if needBuild:
						built = addAtom(a.name, a.element, old_res,
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

def library(cbase, linkers, fragments, rotations=True):
	mol = cbase.atoms[0].molecule
	explore = product(range(len(linkers)), range(len(fragments)))
	library = {}
	for i, j in explore:
		#TODO: Create custom molecule_from_atoms()
		new = SplitMolecule.split.molecule_from_atoms(mol, cbase.atoms)
		target = new.atoms[-1]
		linker = insert(linkers[i], target=target)
		linker_anchor = [ a for a in linker if a.anchor in (4,6,8)][0]
		insert(fragments[j], target=linker_anchor)
		if rotations: 
			fragment_anchor = [ a for a in linker_anchor.neighbors if a.element.number != 1
				and a not in linker ]
			target_neighbor = [ a for a in target.neighbors if a.element.number != 1
				and a not in linker ][0]
			bonds = utils.box.sequential_bonds([target]+linker[:]+fragment_anchor,target_neighbor)
			
			bondrots = []
			for b in bonds:
				br = chimera.BondRot(b)
				br.myanchor = utils.box.find_nearest(new.atoms[0], b.atoms)
				bondrots.append(br)
			library[i,j] = [new, bondrots]
		else:
			library[i,j] = [new]
	return library