import os
import chimera
from chimera import runCommand as rc #cmd line features


## CONFIGURATION ##########################################
# list of atoms involved (defined by user in midas format)
bond_atoms = ":0OD@C11@C12" #torsion angle
bond_atoms2 = ":0OD@C12@C13" #torsion angle
clash_atoms = ":0OD:UNK:HAXSOP" #center of the clashable area
hbond_atoms = ":0OD:UNK:HAXSOP" #normally, the ligand
#torsion step
step = 50
#working directory (create a /tmp dir inside if needed)
wd = "C:/x/"
if not os.path.isdir(wd + "tmp/"):
	os.makedirs(wd + "tmp/")

############################################################
## DONT TOUCH ##############################################
rc("rotation 1 reverse " + bond_atoms) # activate rotation
rc("rotation 2 reverse " + bond_atoms2) # activate rotation
rc("sel " + clash_atoms + " za < 5") #sel everything within 8A
rc("namesel clashable")
rc("sel " + hbond_atoms) #sel everything within 8A
rc("namesel hbonds")

# Results file
output = open(wd + "data.csv", 'w')
output.write("Degree 1_Degree 2;Contacts;H bonds\n")
# Now, loop
for i in xrange(1,361,step): # edit the step, if needed
	rc("rotation 1 " + str(step)) #rotate 'step' degrees
	for j in xrange(1,361,step):
		rc("rotation 2 " + str(step)) #rotate 'step' degrees
		#calculate clashes
		fn_clash = wd + "tmp/degree__" + str(i) + "_" + str(j) + ".clash.txt"
		rc("findclash clashable saveFile " + fn_clash)
		#calculate H bonds
		rc("sel " + hbond_atoms) #sel everything within 8A
		fn_hbonds = wd + "tmp/degree__" + str(i) + "_" + str(j) + ".hbonds.txt"
		rc("hbonds selRestrict any saveFile " + fn_hbonds) #find clashes

		#parse clashes
		f1 = open(fn_clash, 'r')
		contacts = f1.readlines()
		num_of_contacts = contacts[6].split(" ", 1)
		f1.close()
		#parse H bonds; careful! the source needs to be modified for this to work
		#since it will look for a summary line not present in original sc
		f2 = open(fn_hbonds, 'r')
		hbonds = f2.readlines()
		num_of_hbonds = hbonds[0].split(" ", 1)
		f2.close()

		output.write(str(i) + "_" + str(j) + ";" + str(num_of_contacts[0]) + ";" + str(num_of_hbonds[0]) + "\n")
		

output.close()
# clean tmp
[ os.remove(wd + "tmp/" + f) for f in os.listdir(wd + "tmp/") if f.endswith(".txt") ]
os.rmdir(wd + "tmp/")