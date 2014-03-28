# Hyde project
Given an anchored ligand, explores conformational space of randomly generated substituents and optimizes selected residues in target protein with a genetic algorithm.

# Roadmap
## New feats
* ~~Mutamers instead of rotamers~~
* Choose a random H of insertion
* Switch genes on or off on user request
* Improve mutation and crossover functions
	- Allow building blocks evolution (?)

## Performance optimization
* Reduce search scope of clashes and H bonds detection (within `n` angstroms of ligand)
	- Maybe something like ```n = sum(bond.length in ligand)```
* ~~Render all possible molecules at the beginning, rather than every iteration~~
* Review insertMol() looping. It may be redundant / too recursive
* Add multithreading support
* Benchmark the weights

## Known issues
* Reported ~~clashes and~~ H bonds values < actual values!

## Coding
* Clean the code
* Use hierarchical modules
* Convert to OOP