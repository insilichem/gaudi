# Hyde project
Given an anchored ligand, explores conformational space of randomly generated substituents and optimizes selected residues in target protein with a genetic algorithm.

## Usage
In Chimera, define three selections with `namesel`:

* `base`: Single atom. Terminal end of the static part of the ligand

* `anchor`: Single atom. Together with `base`, it will define the static part. It's also the insertion point for the linker.

* `mutable`: Residues to be swapped and mutated.

Then, in command line, type `runscript /abs/path/to/hyde7.py -p <size of initial pop> -g <num of generations>`.

Check the roadmap in tab "Issues".