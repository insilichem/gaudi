#!/usr/bin/python

# MOF3D
# Multi-Objective Force-Field-Free Docking
# A docking module for UCSF Chimera
# Jaime RGP <https://bitbucket.org/jrgp> @ UAB, 2014

import chimera
import os

def atoms_between(atom1, atom2):
	''' Finds all connected atoms between two given atoms '''
	chain1 = [atom1]
	chain2 = [atom2]
	i = 0
	while i < len(chain1):
		a1 = chain1[i]
		if atom2 not in a1.neighbors:
			chain1.extend([ a for a in a1.neighbors if a not in chain1 ])
		i += 1
	j = 0
	while j < len(chain2):
		a2 = chain2[j]
		if atom1 not in a2.neighbors:
			chain2.extend([ a for a in a2.neighbors if a not in chain2 ])
		j += 1
	
	return set(chain1) & set(chain2)

def atoms_by_serial(*serials, **kw):
	if not kw['atoms']:
		kw['atoms'] = [ a for m in chimera.openModels.list() for a in m.atoms ]
	return [a for a in kw['atoms'] if a.serialNumber in serials]

def files_in(path, ext=None):
	if ext:
		return [ os.path.join(path,fn) for fn in next(os.walk(path))[2] if fn.endswith('.'+ext) ]
	return [ os.path.join(path,fn) for fn in next(os.walk(path))[2] ]

def find_nearest(anchor, atoms):
	''' Returns closer atom from `atoms` to `anchor` '''
	try: return next(a for a in atoms if a is anchor)
	except StopIteration: #oops, didn't find it...
		return min(atoms, key= lambda a:len(atoms_between(anchor, a)))

def highest_atom_indices(r):
	''' Returns a dictionary with highest atom indices in given residue 
		Key: value -> element.name: highest index in residue
	'''
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

def pseudobond_to_bond(molecule, remove=False):
	''' Transforms every pseudobond in `molecule` to a covalent bond '''
	pbgroup = chimera.misc.getPseudoBondGroup(
				"coordination complexes of %s (%s)" % 
				(molecule.name, molecule))#, associateWith=[molecule])
	if pbgroup.pseudoBonds:
		for pb in pbgroup.pseudoBonds:
			chimera.molEdit.addBond(*pb.atoms)
			if remove: pbgroup.deletePseudoBond(pb)
		pbm = molecule.pseudoBondMgr()
		pbm.deletePseudoBondGroup(pbgroup)

def write_individuals(inds, outpath, name, evalfn):
	from WriteMol2 import writeMol2
	if not os.path.isdir(outpath):
		os.makedirs(outpath)
	results = []
	for i, ind in enumerate(inds):
		ligand = evalfn(ind, close=False)
		fullname = '{}{}__{:02d}.mol2'.format(outpath, name, i+1)
		writeMol2([ligand], fullname)
		chimera.openModels.remove([ligand])
		results.append([fullname, ind.fitness])
	return results

def sequential_bonds(atoms,s):
	''' Returns bonds in `atoms` in sequential order, beginning at atom `s` '''
	if s not in atoms: atoms.append(s)
	bonds = list(set([ b for a in atoms for b in a.bonds if set(b.atoms)<=set(atoms) ]))
	nbonds = []
	while bonds:
		b = bonds.pop(0)
		if s in b.atoms:
			nbonds.append(b)
			s = b.otherAtom(s)
		else: 
			bonds.append(b)
	return nbonds
