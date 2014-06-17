GAUDI Input files full specification
====================================

GAUDI uses YAML-formatted files for both input and output files. YAML is a human-readable serialization format, already implemented in a broad range of languages. Formally, it consists of a number of dictionaries, whose values are dictionaries themselves. However, due to YAML high readibility, it looks just like a typical indented list.

First-level dictionaries will be referred as "sections", while second-level entries will be called "parameters". GAUDI currently accepts six different sections, five of which are always required:

- default
- ga
- protein
- ligand
- rotamers (optional)
- objectives

## default
This section is meant to store general purpose parameters, such as `savepath` and `savename`.

- **str `savepath`**. Absolute path to the desired output folder, where output files will be generated. If a name conflict arises, new files will overwrite existing ones, so be careful not to choose a directory that contains important information.
- **str `savename`**. A tag prefix that will identify the output files.
- **int `precision`**. Number of decimals that will be printed to stdout progress reports. It has no effect on actual calculations.
- **int `results`**. Number of results to save and write. It has no effect if `pareto` is enabled.
- **bool `pareto`**. If set to `1`, it will return the whole Pareto front, lexicographically sorted. If set to `0`, it will return a list of the best individuals that existed during the evolution, according to its lexicographical scoring. Recommended option: `1`

## ga
- **int `pop`**. Number of individuals in initial population.
- **int `gens`**. Number of generations that will be simulated.
- **float `mu`**.
- **float `lambda_`**.
- **int `mut_eta`**.
- **int `mut_pb`**.
- **int `mut_indpb`**.
- **int `cx_eta`**.
- **int `cx_pb`**.

## protein
- **str `path`**. Absolute path to the protein mol2 file.
- **int `origin`**. Serial number of the protein atom that will serve as center for the exploration sphere.
- **float `radius`**. Radius of the exploration sphere, in angstroms. Please note that this sphere only constrains the position of the `donor` atom of the ligand. Because of that, if the donor atom is located near the terminal ends of the ligand, the effective sphere can be greater that thought.

## ligand
- **str `path`**. Absolute path to either:
    - *A single mol2 file*. Classic docking essay.
    - *A directory that contains mol2 files*. Competitive multiligand docking essay. 
      All the mol2 files found in that path will try to find its best pose in the protein.
    - *A directory that contains subdirectories with mol2 files*. Blocks multiligand essay. 
      GAUDI will build all possible ligands with the molecular units provided in each subdirectory. The subdirectories can be freely named, but do take into account their alphabetic sorting, since GAUDI will respect that ordering for building the ligands. 
      I.e., if you provide a path that contains folders `01/`, `02/`, and `03/`, GAUDI will take a random mol2 from `01/` as the starting point, then append a random mol2 from `02/` and, finally, append a random mol2 from `03/`.
  
    For every mol2 file, an homonimous `attr` file should be provided (`ligand.mol2` should be accompanied by `ligand.attr`). This file must define two main ligand atoms in YAML format, such as:

        atoms:
          acceptor: serialNumber of the ligand atom that will bond to new molecular blocks.
          donor: serialNumber of the ligand atom that will bond to the acceptor of another ligand. It also works as an anchor for the spatial exploration in free docking mode.
- **str `type`**. *deprecated*
- **float `flexibility`**. Number of degrees each rotatable bond can torsion either way. For example, if set to `20`, the bonds will be able to torsion 10 degrees clockwise and 10 degrees counterclockwise. Set it to `360` to allow full free rotation.
- **str `assess`**. *experimental* Absolute path to a ligand pose that will act as RMSD benchmark.

## rotamers (optional)
- **list int `residues`**. List of the positions of the residues that will be able to rotate.
- **str `library`**. Allows to choose between Dunbrack and Dynameomics.
- **int `top`**. Limit eligible rotamers to top n results, based on probability.
- **list of str `mutate`**. If set, a list of three-letter aminoacid codes the residues will be able to mutate into. It can be set to `all` instead of the full list.


## objectives
This section is somehow special, because it features a new level between the parameters and the section itself. It uses dashes for every objective, turning its contents in a list of dictionaries. As a result, the indentation must be strictly respected.

- **str `name`**. Custom name to identify the objective in further analysis. Spaces are not recommendend, but can be used provided you quote the string.
- **str `type`**. Choose between `hbonds`, `contacts`, `distance`, `solvation` and `rmsd`. More objectives will be implemented in newer versions. Each type can have its own parameters.
- **float `weight`**. If `n > 0`, the objective will be maximized. If `n < 0`, it will be minimized. It can accept any float, but it works fine with `1.0` and `-1.0`, respectively. 

Objective-specific parameters are discussed below.

### hbonds
- **list of int `targets`**.
- **float `distance_tolerance`**.
- **float `angle_tolerance`**.

### contacts
- **str `which`**.
- **str `test`**.
- **float `threshold`**.
- **float `threshold_h`**.
- **float `threshold_c`**.
- **float `cutoff`**.


### distance
- **int `target`**.
- **list of int `probes`**.
- **float `threshold`**.
- **float `threshold2`**.

### solvation
- **str `which`**.

### rmsd
- **str `path`**.
