## Create attr file from mol2
import sys

# Open and parse mol2 file
mol2 = sys.argv[1]
attr = sys.argv[2]


out = mol2.replace('.mol2', '.attr')
f = open(mol2, 'r')

atoms = []
parsing = False
for line in f.readlines():
	if line.startswith('@<TRIPOS>ATOM'):
		parsing = True
		continue
	elif line.startswith('@<TRIPOS>BOND'):
		break

	if not parsing:
		continue
	
	words = line.split()
	atoms.append(words[1])
f.close()

f = open(out,'w')

f.write("""# Attr file generated with mol2_to_attr.py
attribute: {0}
match mode: any
recipient: atoms
# list of atoms
""".format(attr))

for a in atoms:
	f.write("\t@{0}\t0\n".format(a))

f.close()