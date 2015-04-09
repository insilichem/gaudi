# GAUDI Project
GAUDI (Genetic Algorithm for Unified Docking Inference) provides a novel method in docking essays by combining several chemostructural criteria, along with energy optimization.

## Features
**True multi-objective optimization**

* Feel free to optimize H bonds, hydrophobic interactions, desolvation effects, distances between given sets of atoms, rotamers and more.

**Force-field-less approach**

* Metal complexes are more than welcome.

**Unprecedented customizability**

* Each objective can be switched off upon demand.

* The weights can be auto-optimized for each essay.

## Installation and dependencies
To download a copy of Gaudi, just clone it to your desired directory using `git clone https://jrgp@bitbucket.org/jrgp/gaudi.git`. It will create a new subdirectory called `gaudi`.

**Dependencies**

* [UCSF Chimera](http://www.cgl.ucsf.edu/chimera/download.html).

* [PyYaml]() package

* [repoze.lru]() package

**Install Chimera**

1. Download the [latest stable copy](http://www.cgl.ucsf.edu/chimera/download.html) and install it with `chmod +x chimera-*.bin && ./chimera-*.bin`. A text-based wizard will guide you. 

2. Locate the main binary. The default location is `~/.local/UCSF-Chimera<arch>-<version>/bin/chimera`, so adapt it for your version. For example, the binay of a 64-bit Chimera v1.10.1 will be at `~/.local/UCSF-Chimera64-1.10.1/bin/chimera`

3. Create a bash alias. Open `~/.bashrc` with your favourite editor and add this line at the end of the file `alias chimera=~/.local/UCSF-Chimera64-1.10.1/bin/chimera`. Save it.

4. Apply the changes with `source ~/.bashrc`. This way, you can run Chimera from the console by typing `chimera`.

**Install Python packages**

UCSF Chimera uses its own copy of Python 2.7, so the packages must be installed with that one, and not the system's. For each package, follow these instructions:

1. Download the source of the package, usually as a *.tar.gz file, and extract it to a new directory with `tar -zxvf <thepackage>.tar.gz`. Once completed, `cd` into that directory.

2. Instead of using the regular `python setup.py install`, you will install it with Chimera's Python: `chimera --nogui --script "setup.py install"`.

3. The bundled Python will copy the contents to the needed folders and precompile the files for you. And that's it!

## Running a GAUDI job
You only have to run `base.py <inputfile>.gaudi` with Chimera's Python. You have already done this to install the packages!

1. `cd /path/to/gaudi/`

2. `chimera --nogui --script "base.py /path/to/input/file.gaudi"`