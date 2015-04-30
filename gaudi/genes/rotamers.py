#!/usr/bin/python

##############
# GAUDIasm: Genetic Algorithms for Universal
# Design Inference and Atomic Scale Modeling
# Author: Jaime Rodriguez-Guerra Pedregal
# Email: jaime.rogue@gmail.com
# Web: https://bitbucket.org/jrgp/gaudi
##############

# Python
import random
from collections import OrderedDict
# Chimera
from Rotamers import getRotamers, useRotamer, NoResidueRotamersError
import SwapRes
# External dependencies
from repoze.lru import LRUCache
import deap.tools
# GAUDI
from gaudi.genes import GeneProvider
from gaudi.parse import parse_rawstring

class Rotamers(GeneProvider):

	def __init__(self, parent=None, name=None, cache=None,
				residues=None, library='Dunbrack', mutations=[],
				**kwargs):
		self._kwargs = kwargs
		self.name = name
		self.parent = parent
		self._cache = cache
		self._residues = residues
		self.library = library
		self.mutations = mutations
		self.allele = []
		# set (or retrieve) caches
		if self.name not in self._cache:
			self._cache[self.name] = {	'residues': OrderedDict(),
										'rotamers': LRUCache(300) }
		self.residues = self._cache[self.name]['residues']
		self.rotamers = self._cache[self.name]['rotamers']
		
		# find requested residues
		self._residues_rawstring = tuple(parse_rawstring(r) for r in residues)
		for molecule, resid in self._residues_rawstring:
			try:
				res = next(r for r in self.parent.genes[molecule].compound.mol.residues
								if r.id.position == resid)
			except (KeyError, StopIteration): # molecule or residue not found
				raise
			else: #residue was found!
				self.residues[(molecule, resid)] = res
				self.allele.append(
					(	random.choice(self.mutations+[res.type]),
						random.random()
					)
				)

	def __deepcopy__(self, memo):
		new = self.__class__(self.parent, self.name, self._cache,
							self._residues, self.library, self.mutations,
							**self._kwargs)
		new.__dict__.update((k,v) for k,v in self.__dict__.items())
		new.allele = self.allele[:]
		return new 

	def express(self):
		for (mol, pos), (restype, i) in zip(self.residues, self.allele):
			try:
				rot = self.get_rotamers(mol, pos, restype)
			except NoResidueRotamersError: # ALA, GLY...
				SwapRes.swap(self.residues[(mol,pos)], restype)
			else:
				useRotamer(self.residues[(mol,pos)], [rot[int(i*len(rot))]])	
			
	def unexpress(self):
		pass

	def mate(self, mate):
		self.allele, mate.allele = deap.tools.cxTwoPoint(self.allele, mate.allele)

	def mutate(self, indpb):
		self.allele = []
		for molecule, resid in self._residues_rawstring:
			try:
				res = next(r for r in self.parent.genes[molecule].compound.mol.residues
								if r.id.position == resid)
			except (KeyError, StopIteration): # molecule or residue not found
				raise
			else: #residue was found!
				self.residues[(molecule, resid)] = res
				self.allele.append(
					(	random.choice(self.mutations+[res.type]),
						random.random()
					)
				)

	def write(self, path, name):
		pass
		# rotamerline = '{}.{} {} {} {}\n'.format(res.id.position, res.id.chainId,
		# 			self.library.title(), res.type, ' '.join(map(str,rot.chis)))

	###
	def get_rotamers(self, mol, pos, restype):
		rotamers = self.rotamers.get((mol, pos, restype))
		if not rotamers:
			try:
				rotamers = getRotamers(self.residues[(mol,pos)], resType=restype, 
										lib=self.library.title())[1]
			except NoResidueRotamersError: # ALA, GLY... has no rotamers
				raise
			except KeyError:
				raise
			else:
				self.rotamers.put((mol,pos,restype), rotamers)
		return rotamers

def enable(**kwargs):
	return Rotamers(**kwargs)