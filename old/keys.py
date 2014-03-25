import chimera

def countClashes(atoms=chimera.selection.currentAtoms()):
	#calculate clashes
	import DetectClash as dc
	clashes = dc.detectClash(atoms,intraRes=True)
	num_of_clashes = 0
	if clashes:
		for clashList in clashes.values():
			num_of_clashes += len(clashList)

	return clashes


def allKeys(d):
	for k, v in d.iteritems():
		if isinstance(v, dict):
			for k_ in allKeys(v): 
				yield k_
			yield k
		else:
			yield k



dic = { 0: { 1 : { 3: 4}}, 2: 2 }
print dic
clashes = countClashes()
#c = allKeys(clashes)
for c in allKeys(clashes): print c