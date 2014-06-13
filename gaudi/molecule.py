#!/usr/bin/python

# gaudi
# Genetic Algorithm for Unified Docking Inference
# A docking module for UCSF Chimera
# Jaime RGP <https://bitbucket.org/jrgp> @ UAB, 2014

import chimera, BuildStructure
from chimera import UserError
from chimera.molEdit import addAtom, addBond
import os
import itertools
import move, utils

class Library(object):

	def __init__(self, path, origin=None, covalent=False, flexible=True):
		self.path = path
		self.origin = origin if origin else chimera.Point(0,0,0)
		self.covalent = covalent
		self.flexible = flexible
		self.compounds = {}
		self.compile_catalogue()
	
	def __getitem__(self, index):
		return self.get(index)
	
	def get(self, index):
		try:
			return self.compounds[index]
		except KeyError:
			self.compounds[index] = self.build(index, self.origin)
			return self.compounds[index]
	
	def compile_catalogue(self):
		if os.path.isdir(self.path):
			folders = sorted([ os.path.join(self.path, d) for d in os.listdir(self.path)
					if os.path.isdir(os.path.join(self.path,d)) and not d.startswith('.')
						and not d.startswith('_')])
			if folders:
				self.catalog = list(itertools.product(*[utils.box.files_in(f, ext='mol2')
														 for f in folders]))
			else:
				self.catalog = [ (f,) for f in utils.box.files_in(self.path, ext='mol2') ]
		elif os.path.isfile(self.path) and self.path.endswith('.mol2'):
			self.catalog = [(self.path,)]

	def build(self, index, where=None):
		if self.covalent:
			base = Compound()
			base.donor = base.add_dummy_atom(where.neighbors[0], serial=1)
			base.acceptor = base.add_dummy_atom(where, bonded_to=base.donor, serial=2)
			base.append(Compound(molecule=index[0]))
		else:
			base = Compound(molecule=index[0])
			base.place(where)

		for molpath in index[1:]:
			base.append(Compound(molecule=molpath))
		
		if self.flexible:
			base.update_rotatable_bonds()
		return base

class Compound(object):
	### Initializers
	def __init__(self, molecule=None, origin=None, **kwargs):
		if isinstance(molecule, chimera.Molecule):
			self.mol = molecule
		elif not molecule or molecule=='dummy':
			self.mol = _dummy_mol('dummy')
		else:
			self.mol = chimera.openModels.open(molecule)[0]
			chimera.openModels.remove([self.mol])
		
		self.mol.gaudi = self
		self.rotatable_bonds = []
		self.parse_attr()
		self.origin = origin
		for k,v in kwargs.items():
			self.__dict__[k] = v
		if not hasattr(self, 'nonrotatable'):
			self.nonrotatable = []

	def parse_attr(self):
		try:
			attr = utils.parse.Settings(self.mol.openedAs[0][:-4]+'attr', asDict=True).parsed
		except (AttributeError, IOError):
			print "No attr file found"
		else:
			flat_attr = {}
			if 'atoms' in attr:
				for k,v in attr['atoms'].items():
					flat_attr[k] = next(a for a in self.mol.atoms if a.serialNumber==v)
			if 'angles' in attr:
				flat_attr['angles'] = {}
				for k,v in attr['angles'].items():
					flat_attr['angles'][tuple(k)] = v
			if 'bonds' in attr:
				for k,v in attr['bonds'].items():
					if k == 'nonrotatable':
						if isinstance(v, list):
							flat_attr[k] = utils.box.atoms_by_serial(*v, atoms=self.mol.atoms)
						elif v.strip().lower() in ('all', 'yes'):
							flat_attr[k] = self.mol.atoms

			self.__dict__.update(flat_attr)

	def update_attr(self, d):
		for k,v in self.__dict__.items():
			if k in ('mol', 'angles'):
				continue
			if isinstance(v, list):
				setattr(self, k, [d[v_] if v_ in d else v_ for v_ in v])
			else:
				setattr(self, k, d[v] if v in d else v)
	
	def get_rotatable_bonds(self):
		existing_bondrots = self.rotatable_bonds
		existing_bondrots_bonds = []
		for br in existing_bondrots:
			existing_bondrots_bonds.append(br.bond)
			yield br

		bonds = set(b for a in self.mol.atoms for b in a.bonds if not a.element.isMetal)
		bonds = sorted(bonds, key=lambda b: min(y.serialNumber for y in b.atoms))

		for b in bonds:
			if b in existing_bondrots_bonds:
				continue
			a = b.atoms[0]
			if  a not in self.nonrotatable and \
				a.idatmType in ('C3', 'N3', 'C2', 'N2') and \
				(a.numBonds > 1 and b.otherAtom(a).numBonds > 1) or \
				a.name == 'DUM' or b.otherAtom(a).name == 'DUM' :
					try:
						br = chimera.BondRot(b)
					except (chimera.error, ValueError), v:
						if "cycle" in str(v):
							continue #discard bonds in cycles!
						else:
							raise
					else:
						br.rotanchor = utils.box.find_nearest(self.donor, b.atoms)
						yield br

	def update_rotatable_bonds(self):
		self.rotatable_bonds = list(self.get_rotatable_bonds())

	def destroy(self):
		chimera.openModels.close([self.mol])
		del self

	####
	# Building and moving functions
	####
	def add_dummy_atom(self, where, name='dum', element=None, residue=None,
						bonded_to=None, serial=None):
		if isinstance(where, chimera.Atom):
			element = where.element if not element else element
			where = where.coord()
		else:
			element = chimera.Element('C')
		
		residue = self.mol.residues[-1] if not residue else residue
		return addAtom(name, element, residue, where, serial, bonded_to)

	def append(self, molecule):
		self.attach(molecule, self.acceptor, molecule.donor)

	def prepend(self, molecule):
		self.attach(molecule, self.donor, molecule.donor)

	def attach(self, molecule, acceptor, donor):
		if acceptor not in self.mol.atoms:
			raise UserError('Specified atom is not part of molecule.')
		 
		molecule.place_for_bonding(acceptor)
		if not donor:
			donor = molecule.donor
		updated_atoms = self.join(molecule, acceptor, donor)

		# Update Gaudi ATTR
		self.update_attr(updated_atoms) #update existant atoms
		self.acceptor = updated_atoms[molecule.acceptor] # change acceptor
		self.axis_end = updated_atoms[molecule.axis_end]
		if hasattr(molecule, 'nonrotatable'):
			nonrot_atoms = [updated_atoms[a] for a in molecule.nonrotatable]
			if hasattr(self, 'nonrotatable'):
				self.nonrotatable.extend(nonrot_atoms)
			else:
				self.nonrotatable = nonrot_atoms

		molecule.destroy()

	def join(self, molecule, acceptor, donor, newres=False):
		target = acceptor
		sprouts = [ donor ]
		res = _dummy_res() if newres else target.residue
		res_atoms = res.atoms

		i = max(a.serialNumber for a in self.mol.atoms)
		index = utils.box.highest_atom_indices(self.mol)
		for a in molecule.mol.atoms:
			try:
				index[a.element.name] += 1
			except KeyError:
				index[a.element.name] = 1
			finally:
				a.name = a.element.name + str(index[a.element.name])
		built_atoms = {}		
		while sprouts:
			sprout = sprouts.pop(0) # get first atom
			if sprout.name in res.atomsMap:
				target = res.atomsMap[sprout.name][-1]
			else:
				i += 1
				built_atoms[sprout] = target = addAtom(sprout.name, 
								sprout.element, res, sprout.coord(),
								bondedTo=target, serialNumber=i)

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
					built_atoms[a] = built = addAtom(a.name, a.element, res, a.coord(), 
									bondedTo=target, serialNumber=i)
					# if a has more than one neighbor:
					if len(a.neighbors) > 1:
						sprouts.append(a) # this new atom can be a new sprout
				if built not in target.bondsMap and built not in res_atoms: #link!
					addBond(target, built)

		return built_atoms

	# def orientate(self, target, **kwargs):
	# 	if kwargs:
	# 		angles = kwargs
	# 	elif hasattr(self, 'angles'):
	# 		angles = self.angles
	# 	else:
	# 		raise UserError('Please specify geometric constraints for orientation.')

	# 	for atoms, alpha in angles.items():
	# 		parsed_atoms = []
	# 		post_atoms = []
	# 		for a in atoms: #parse atoms
	# 			if isinstance(a, chimera.Atom):
	# 				parsed_atoms.append(a)
	# 			a = a.strip()
	# 			if a == 'post':
	# 				post_atoms = self.donor.neighbors
	# 			elif a == 'anchor':
	# 				parsed_atoms.append(self.donor)
	# 			elif a == 'target':
	# 				parsed_atoms.append(target)
	# 			elif a == 'pre':
	# 				parsed_atoms.append(next(a for a in target.neighbors if a.numBonds>1))
	# 			elif a == 'axis':
	# 				parsed_atoms.extend([self.axis_start, self.axis_end])
	# 			elif isinstance(a, str):
	# 				raise chimera.UserError("Atom descriptor {} not supported.".format(a))
			
	# 		for pa in post_atoms:
	# 			move.rotate(self.mol, parsed_atoms+[pa], alpha)

	# def orientate(self, anchor, target):
	# 	target_geometry = target

	def place(self, where, anchor=None):
		if isinstance(where, chimera.Atom):
			where = where.coord()
		if not anchor:
			anchor = self.donor
		move.translate(self.mol, anchor, where)

	def place_for_bonding(self, target, anchor=None):
		if not isinstance(target, chimera.Atom):
			raise UserError('Target must be a chimera.Atom object.')
		if not anchor:
			anchor = self.donor

		
		def new_atom_position(atom, newelement):
			geometry = chimera.idatm.typeInfo[atom.idatmType].geometry
			bond_length = chimera.Element.bondLength(atom.element, newelement)
			neighbors_crd = [a.coord() for a in atom.neighbors]
			return chimera.bondGeom.bondPositions(atom.coord(), geometry, bond_length,
					neighbors_crd)[0]
		# Get target position
		target_pos = new_atom_position(target, anchor.element)
		# Place it
		self.place(target_pos)
		# Fix orientation
		anchor_pos = new_atom_position(anchor, target.element)
		move.rotate(self.mol, [target.coord(), anchor.coord(), anchor_pos], 0.0)
		# # Hydrogen is no longer needed
		# target.molecule.deleteAtom(hydrogen)


def _add_hydrogen(atom):
	h = None
	if atom.element.number == 1:
		return atom.neighbors[0], atom
	try:
		geometry = chimera.idatm.typeInfo[atom.idatmType].geometry
	except KeyError:
		geometry = 3
		print "Warning, arbitrary geometry with {}, {}".format(atom, geometry) 
	
	try:
		atom, h = BuildStructure.changeAtom(atom, atom.element, 
								geometry, atom.numBonds+1)[:2]
	except ValueError:
		try: 
			h = next(a for a in atom.neighbors if a.numBonds == 1)
		except StopIteration:
			h = next(a for a in atom.molecule.atoms if a.numBonds == 1)
			print "WARNING! Atom {} was chosen randomly among molecule terminal atoms".format(h)
		else:
			print "WARNING! Atom {}, with geometry {} and {} bonds was chosen randomly among neighbor terminal atoms".format(
												atom, geometry, atom.numBonds)		
	return atom, h

def _dummy_mol(name):
	m = chimera.Molecule()
	m.name = name
	r = m.newResidue(name, 'het', 1, ' ')
	r.isHet = True
	return m

		# # Create a new hydrogen
		# target_atom, hydrogen = _add_hydrogen(target_atom)
		# # Adjust bond length to fit new bond
		# dv = hydrogen.coord() - target_atom.coord()
		# dv.length = chimera.Element.bondLength(self.donor.element, target_atom.element)
		# diff = dv - (hydrogen.coord() - target_atom.coord())
		# hydrogen.setCoord(hydrogen.coord() + diff)