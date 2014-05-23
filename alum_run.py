import chimera
import glob, sys
pdbs = glob.glob(sys.argv[1]+"*.pdb")

for pdb in pdbs:
	mol = chimera.openModels.open(pdb)
	chimera.runCommand('del element.H')
	chimera.runCommand('runscript /home/jr/x/hyde/alum.py -p 200 -g 10 -t 2.0')
	chimera.openModels.close(mol)
