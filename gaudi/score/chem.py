#!/usr/bin/python

# gaudi
# Genetic Algorithm for Unified Docking Inference
# A docking module for UCSF Chimera
# Jaime RGP <https://bitbucket.org/jrgp> @ UAB, 2014

import chimera, FindHBond, DetectClash, Measure, MoleculeSurface, ChemGroup as cg

ALIPH = ['C3', [[cg.C, [cg.C, [cg.C , [cg.R, cg.R, cg.R, cg.R]], \
		cg.R, cg.R] ], cg.R, cg.R, cg.R] ], [ 1, 1, 1, 1, 1, 0, 0] 
PI = 3.14159265359

# Callables
def hbonds(models=None, sel=None, selRestrict=True, cache=False, test=None, 
			dist_slop=0.4, angle_slop=20.0):
	## Calculates H bonds in `model` given `sel` atoms
	if not models:
		models = chimera.openModels.list(modelTypes=[chimera.Molecule])
	if not sel and selRestrict:
		sel = chimera.selection.currentAtoms()
	#calculate H bonds
	hbonds = FindHBond.findHBonds(models, cacheDA=cache, donors=test, acceptors=test,
									distSlop=dist_slop, angleSlop=angle_slop)
	if selRestrict:
		try:
			hbonds = FindHBond.base._filterBySel(hbonds, sel, 'any')
		except AttributeError: #renamed in daily builds!
			hbonds = FindHBond.base.filterHBondsBySel(hbonds, sel, 'any')
	
	return hbonds

def clashes(atoms, test='others', clashThreshold=0.6,
		hbondAllowance=0.4, assumedMaxVdw=2.1,
		bondSeparation=4, intraRes=True, interSubmodel=True,
		parse=False, parse_threshold=0.4):
	#calculate clashes

	clashes = DetectClash.detectClash(atoms, test=test, intraRes=intraRes,
		interSubmodel=interSubmodel, clashThreshold=clashThreshold,
		hbondAllowance=hbondAllowance, assumedMaxVdw=assumedMaxVdw,
		bondSeparation=bondSeparation)
	num_of_clashes = 0
	if clashes:
		num_of_clashes = sum(len(cl) for cl in clashes.values())/2
		if parse:
			pos, neg = _parse_clashes_c(clashes, parse_threshold)
			return clashes, num_of_clashes/2, pos, neg
	#else
	return clashes, num_of_clashes, [], []

def solvation(atoms):
	xyzr_data = Measure.measure.atom_xyzr(atoms)
	surfaces = MoleculeSurface.xyzr_surface_geometry(xyzr_data)
	# return atoms, ses, sas
	return atoms, surfaces[3][:,0], surfaces[3][:,1]

def draw_interactions(interactions, startCol='FF0000', endCol='FFFF00',
		key=None, name="Custom pseudobonds"):
	if not len(interactions):
		return
	pb = chimera.misc.getPseudoBondGroup(name)
	color = _hex_to_rgb(startCol)+[1.0]
	if key != None:	
		max_ = max(abs(_[3]) for _ in interactions)
	for i in interactions:
		npb = pb.newPseudoBond(i[0], i[1])
		if key != None:
			intensity = (max_ - abs(i[key]))/(max_)
			opacity = 1-0.7*intensity
			if startCol != endCol:
				color = _linear_color(intensity, startCol, endCol)+[opacity]
			else:
				color = _hex_to_rgb(startCol)+[opacity]
		npb.color = chimera.MaterialColor(*color)

## Internal use
def _hex_to_rgb(hexa):
	return [ int(hexa[i:i+2], 16) for i in range(0,6,2) ]
def _linear_color(value, start, end):
	color = []
	c = 0
	for s, e in zip(_hex_to_rgb(start),_hex_to_rgb(end)):
		s, e = sorted([s, e])
		c = s + value*(e-s)
		color.append(c/255.)
	return color

def _lennard_jones(a1, a2):
	dist = a1.xformCoord().distance(a2.xformCoord())
	zero = 0.98*(a1.radius + a2.radius)
	x = zero/dist
	return (x**12 - 2*x**6)

def _parse_clashes(clashes, mol, threshold):
	m = chimera.openModels.list(modelTypes=[chimera.Molecule])
	AROMATIC = set(a for g in cg.findGroup("aromatic ring",m) for a in g)
	ALIPHATIC = set(a for g in cg.findGroup(ALIPH, m) \
					for a in g if a not in AROMATIC)

	positive, negative = [], []
	for a1, c in clashes.items():
		for a2, dist in c.items():
			if dist <= threshold and a1.residue != a2.residue:
				if a1 in AROMATIC and a2 in AROMATIC:
					positive.append([a1, a2, dist, _lennard_jones(a1, a2)])
				elif a1 in AROMATIC and a2 in ALIPHATIC:
					positive.append([a1, a2, dist, _lennard_jones(a1, a2)])
				elif a1 in ALIPHATIC and a2 in AROMATIC:
					positive.append([a1, a2, dist, _lennard_jones(a1, a2)])
			elif dist > threshold:
				negative.append([a1, a2, dist, _vdw_vol_overlap(a1,a2)])

	return positive, negative

def _parse_clashes_c(clashes, threshold=0.4):
	mols = chimera.openModels.list(modelTypes=[chimera.Molecule])
	vdwatoms = set(a for m in mols for a in m.atoms if a.element.name in ('C', 'S'))

	positive, negative = [], []
	for a1, c in clashes.items():
		for a2, dist in c.items():
			if dist <= threshold and a1.molecule != a2.molecule:
				if a1 in vdwatoms and a2 in vdwatoms:
					positive.append([a1, a2, dist, _lennard_jones(a1, a2)])
			elif dist > threshold:
				negative.append([a1, a2, dist, _vdw_vol_overlap(a1,a2)])

	return positive, negative

def _vdw_vol_overlap(a1, a2):
	# Adapted from Eran Eyal, Comput Chem 25: 712-724, 2004
	d = a1.xformCoord().distance(a2.xformCoord())
	if not d: 
		return 1000
	h_a, h_b = 0, 0
	if d and d < a1.radius+a2.radius:
		h_a = (a2.radius**2 - (d- a1.radius)**2)/(2*d)
		h_b = (a1.radius**2 - (d- a2.radius)**2)/(2*d)
	v = (PI/3) * (h_a**2) * (3*a1.radius - h_a) + \
		(PI/3) * (h_b**2) * (3*a2.radius - h_b)
	return v