import chimera
from BuildStructure import Smiles as s
def buildSprout(template, anchor):
	from chimera.molEdit import addAtom, addBond
	# molecule will always have a prefix C
	sprouts = template.atomsMap['C1']
	builtSprout = [anchor]
	while sprouts:
		sprout = sprouts.pop()
		new_sprout = addAtom(sprout.name, sprout.element, 
			anchor.residue, sprout.coord(), None, bondedTo=builtSprout[-1])
		builtSprot.append(new_sprout)
		for a, b in sprout.bondsMap.items():
			builtBranch = addAtom(a.name, a.element, new_sprout.residue,
				a.coord(), None, bondedTo=new_sprout)


anchor = chimera.selection.currentAtoms()[0]
template = s.smiles2mol("CCC",resName=anchor.residue.type)