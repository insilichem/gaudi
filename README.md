# Hyde project
Some description here

ROADMAP
=======
## New feats
* ~~Mutamers instead of rotamers~~
	- ~~Add one more gene to individual~~
* Choose a random H of insertion

## Performance
* Render all possible molecules at the beginning, rather than every iteration [GHOST branch]
* Review insertMol() looping. It may be redundant
* Add multithreading support
* Benchmark the weights

## Known issues
* Reported ~~clashes and~~ H bonds values < actual values!


## Coding
* Clean the code
	- Variable names, especially `target` and `anchor`
* OOP conversion
