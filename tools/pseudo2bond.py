## Pseudobonds to covalent bonds
import chimera
from chimera import replyobj
from chimera.molEdit import addBond
from chimera.misc import getPseudoBondGroup 
import sys

def pseudo2bond(mol, undo=False):
	pseudobonds = getPseudoBondGroup("coordination complexes of %s (%s)" % 
			(mol.name, mol), associateWith=[mol]).pseudoBonds
	if pseudobonds:
		if undo:
			for pb in pseudobonds:
				a1, a2 = pb.atoms[:2]
				bond = [ b for b in a1.bonds if b in a2.bonds ]
				mol.deleteBond(bond[0])
			return
		#else	
		for pb in pseudobonds:
			addBond(pb.atoms[0],pb.atoms[1])
	else:
		replyobj.status("No pseudobonds exist in current molecule.")

rm = False
if len(sys.argv)>1 and sys.argv[1] == "rm":
	rm = True

mols = chimera.selection.currentMolecules()
if rm:
	[ pseudo2bond(mol, undo=True) for mol in mols ]
else: 
	[ pseudo2bond(mol) for mol in mols ]		