#!/usr/bin/python

# gaudi
# Genetic Algorithm for Unified Docking Inference
# A docking module for UCSF Chimera
# Jaime RGP <https://bitbucket.org/jrgp> @ UAB, 2014

import chimera, BuildStructure
from chimera.molEdit import addAtom, addBond
from itertools import product
import os
#custom
import move, utils

def library(path, bondto=None, join=False, rotations=False):
	default_bondto, default_join = bondto, join
	if os.path.isdir(path):
		folders = sorted([ os.path.join(path, d) for d in os.listdir(path)
					if os.path.isdir(os.path.join(path,d)) and not d.startswith('.') ])
		explore = product(*[utils.box.files_in(f, ext='mol2') for f in folders])
	elif path.endswith('.mol2'):
		explore = [[path]]

	library, new = {}, []
	for x in explore:
		# Build ligands
		bondto, join = default_bondto, default_join
		for mol in x:
			new = place(mol, target=bondto, join=join, inplace=True)
			bondto, join = new.cfg.atoms['target'], True
		# Add rotable bonds, if requested
		if rotations: 
			bondrots = []
			rotanchor = min(new.atoms, key=lambda a:a.serialNumber)
			for b, br in _rotable_bonds(new, sort=True):
				br.rotanchor = utils.box.find_nearest(rotanchor, b.atoms)
				bondrots.append(br)
			library[id(new)] = new, bondrots
		else:
			library[id(new)] = new, None
		chimera.openModels.remove([new])
	return library

def place(mol2, target=None, join=True, p2b=True, inplace=True, geom=None):

	tmpl = chimera.openModels.open(mol2, type="Mol2", hidden=False, shareXform=True)[0]
	_add_attr(mol2.replace('.mol2','.attr'), tmpl)

	if not geom and hasattr(tmpl, 'cfg') and hasattr(tmpl.cfg, 'angles'):
		geom = tmpl.cfg.angles
	
	if inplace:
		anchor = tmpl.cfg.atoms['anchor']
		if isinstance(target, chimera.Point): # free docking mode
			# move.translate(tmpl, tmpl.bbox()[1].center(), target)
			move.translate(tmpl, anchor, target)
		else: #is atom, covalent mode is on
			axis_start, axis_end = tmpl.cfg.atoms['axis_start'],tmpl.cfg.atoms['axis_end']

			#discard H atom and set actual target
			if target.element.number == 1:
				targetH, target = target, target.neighbors[0]
			else: #add new H based on target atom geometry
				try:
					geometry = chimera.idatm.typeInfo[target.idatmType].geometry
				except KeyError:
					print "Warning, arbitrary geometry with {}, {}, {}".format(
															target, geometry)
					geometry = 3
				try:
					target, targetH = BuildStructure.changeAtom(
										target, target.element, geometry, 
										target.numBonds + 1)[:2]
				except ValueError: #unpacking error; pick any terminal atom available
					print "Warning, arbitrary terminal atom with {}, {}, {}".format(
												target, geometry, target.numBonds)
					targetH = next(a for a in target.neighbors if a.numBonds == 1)
					
			# fix bond length
			dv = targetH.coord() - target.coord()
			dv.length = chimera.Element.bondLength(anchor.element, target.element)
			diff = dv - (targetH.coord() - target.coord())
			targetH.setCoord(targetH.coord() + diff)
			H_coord = targetH.coord()
			target.molecule.deleteAtom(targetH)	
			# align target+anchor
			move.translate(tmpl, anchor, H_coord)
			
			if geom:
				for atoms, alpha in geom.items():
					new_atoms = []
					for i, a in enumerate(atoms): #parse atoms
						if isinstance(a, chimera.Atom):
							continue
						elif a == 'post':
							new_atoms.append(anchor.neighbors[0])
						elif a == 'anchor':
							new_atoms.append(anchor)
						elif a == 'target':
							new_atoms.append(target)
						elif a == 'pre':
							new_atoms.append(next(a for a in target.neighbors if a.element.number!=1))
						elif a == 'axis':
							new_atoms.extend([axis_start, axis_end])
						elif isinstance(a, str):
							raise chimera.UserError("Atom descriptor {} not supported".format(a))
					del geom[atoms]
					geom[tuple(new_atoms)] = alpha
					move.rotate(tmpl, new_atoms, alpha)

	if join: # the only way is to create a copy
		return copy_atoms(tmpl.atoms, bondto=target, join=join, keepattr='cfg', close=True)

	return tmpl

def copy_atoms(atoms, bondto=None, join=False, keepattr=None, close=False):
	from collections import OrderedDict

	tmpl = atoms[0].molecule
	utils.box.pseudobond_to_bond(tmpl, remove=True)
	built_atoms = OrderedDict()

	if bondto and join is True:
		target = bondto
		res = target.residue
		mol = target.molecule
		
		i = max(a.serialNumber for a in res.molecule.atoms)
		index = utils.box.highest_atom_indices(res)
		for a in atoms:
			try:
				index[a.element.name] += 1
				a.name = a.element.name + str(index[a.element.name])
			except KeyError:
				index[a.element.name] = 1
				a.name = a.element.name + str(index[a.element.name])
	else:
		i = 0
		res = _dummy_res('ligand')
		mol = res.molecule
		if join == 'dummy':
			built_atoms[bondto] = target = \
				addAtom('DUM', bondto.element, res,
						bondto.coord(), bondedTo=None, serialNumber=2)
			built_atoms[bondto.neighbors[0]] = \
				addAtom('DUM', bondto.neighbors[0].element, res,
						bondto.neighbors[0].coord(), bondedTo=target, serialNumber=1)
			i = 2
	try:
		sprouts = [ tmpl.cfg.atoms['anchor'] ] # start to grow from seed
	except (AttributeError, KeyError):
		print "Warning. Couldn't locate molecular anchor. Starting from arbitrary choice."
		sprouts = atoms[:1]

	res_atoms = res.atoms
	while sprouts:
		sprout = sprouts.pop(0) # get first atom
		if sprout.name in res.atomsMap:
			target = res.atomsMap[sprout.name][-1]
		else:
			i += 1
			built = addAtom(sprout.name, sprout.element, res, sprout.coord(),
							bondedTo=target, serialNumber=i)
			target = built
			built_atoms[sprout] = built
		for a in sprout.neighbors:
			if a.element.number == 1: continue
			if a.name not in res.atomsMap:
				needBuild = True
			else:
				# atom is already present, but it can be part of a cycle
				# if we get to it it's because another atom is linking it
				needBuild = False
				built = res.atomsMap[a.name][-1]
			if needBuild:
				i +=1
				built = addAtom(a.name, a.element, res,	a.coord(), 
								bondedTo=target, serialNumber=i)
				built_atoms[a] = built
				# if a has more than one neighbor:
				if len(a.neighbors) > 1:
					sprouts.append(a) # this new atom can be a new sprout
			if built not in target.bondsMap and built not in res_atoms: #link!
				addBond(target, built)
	
	# Save requested attributes
	if keepattr and hasattr(tmpl, str(keepattr)):
		setattr(mol, keepattr, getattr(tmpl, keepattr))
		if keepattr == 'cfg':
			delattr(mol.cfg, 'angles') #not needed anymore, and outdated
			for k, v in mol.cfg.atoms.items():
				mol.cfg.atoms[k] = built_atoms[v]
			if join == 'dummy':
				mol.cfg.atoms['anchor'] = built_atoms[bondto]
	elif isinstance(keepattr, chimera.Atom):
		mol.kept = built_atoms[keepattr]

	if close:
		chimera.openModels.close(tmpl)
	return mol

#####################################
def _add_attr(attr_file, mol):
	cfg = utils.parse.Settings(attr_file, asDict=True)
	for k, v in cfg.parsed['atoms'].items():
		cfg.parsed['atoms'][k], = utils.box.atoms_by_serial(v, atoms=mol.atoms)
	setattr(mol, 'cfg', utils.parse.Settings.Param(cfg.parsed))

def _center_of_mass(mol):
	# COM = (1/total_mass)*sum(mass_i*coord_i)
	sum_, center = 0, 0
	for a in mol.atoms:
		sum_ += a.element.mass
		center += a.element.mass * chimera.Vector(*a.coord().data())
	mol.com = sum_/center
	return mol.com


def _dummy_res(name, atom=None):
	m = chimera.Molecule()
	m.name = name
	pos = 1
	while m.findResidue(chimera.MolResId('het', pos)) \
	or m.findResidue(chimera.MolResId('water', pos)):
		pos += 1
	r = m.newResidue(name, 'het', pos, ' ')
	r.isHet = True
	chimera.openModels.add([m])
	return r

def _rotable_bonds(mol, sort=False):
	bonds = set(b for a in mol.atoms for b in a.bonds if not a.element.isMetal)
	if sort:
		bonds = sorted(bonds, key=lambda x: min(y.serialNumber for y in x.atoms))
	for b in bonds:
		if any(	a.idatmType in ('C3', 'N3', 'C2', 'N2') and
				(all([b.otherAtom(a).numBonds > 1, a.numBonds > 1]) or
				any([a.name == 'DUM', b.otherAtom(a).name == 'DUM'])) 
				for a in b.atoms
			):
			try:
				br = chimera.BondRot(b)
			except (chimera.error, ValueError), v:
				if "cycle" in str(v):
					continue #discard bonds in cycles!
				else: 
					raise
			else:
				yield b, br