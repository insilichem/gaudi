# findAnchor

from chimera import UserError
from chimera import selection as sel

def atomsBetween(atom1, atom2):
	chain1 = atom1.neighbors
	chain2 = atom2.neighbors
	i = 0
	while i < len(chain1):
		a1 = chain1[i]
		if atom2 not in a1.neighbors:
			chain1.extend( [ a for a in a1.neighbors if a not in chain1 ])
		i += 1
	j = 0
	while j < len(chain2):
		a2 = chain2[j]
		if atom1 not in a2.neighbors:
			chain2.extend( [ a for a in a2.neighbors if a not in chain2 ])
		j += 1
	chain = set(chain1) & set(chain2)
	ends = set([atom1, atom2])

	return chain|ends


def findNearest(anchor, atoms):
	# return atom with mininum connectivity index
	from chimera import Atom

	nearest = Atom
	minimum = 10000

	for a in atoms:
		distance = len(atomsBetween(anchor, a))
		print a, distance
		if distance < minimum:
			minimum = distance
			nearest = a

	return nearest
		
atom1, atom2 = sel.currentAtoms()
# sel.setCurrent(atomsBetween(atom1, atom2))

anchor = sel.currentAtoms()
bonds = sel.currentBonds()

#if len(anchor) > 1: raise UserError("Only one atom can be selected")
#if not anchor in bonds[0].atoms : raise UserError("Anchor atom must belong to a bond")

#anchorSide = findNearest(anchor[0], bonds[0].atoms)
sel.setCurrent(atomsBetween(atom1, atom2))