# GAUDI Project
GAUDIasm: Genetic Algorithms for Universal Design Inference and Atomic Scale Modeling provides a novel method in design essays by combining several chemostructural criteria, along with energy optimization.

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

* **[UCSF Chimera](http://www.cgl.ucsf.edu/chimera/download.html)**. Main framework.

* **[deap](https://pypi.python.org/pypi/deap/) package**. Genetic Algorithms.

* **[PyYaml](http://pyyaml.org/wiki/PyYAML) package**. Input/output file parsing.

* **[repoze.lru](https://pypi.python.org/pypi/repoze.lru/) package**. A lightweight cache.

**Install Chimera**

1 - Download the [latest stable copy](http://www.cgl.ucsf.edu/chimera/download.html) and install it with:

    chmod +x chimera-*.bin && ./chimera-*.bin

2 - Locate the installation directory. The default location is `~/.local/UCSF-Chimera<arch>-<version>`, so adapt it for your version. For example, the binary of a 64-bit Chimera v1.10.1 will be at`~/.local/UCSF-Chimera64-1.10.1`

3 - Now, create some bash aliases to speed up the work. Open `~/.bashrc` with your favourite editor and add these lines at the end of the file.

    CHIMERADIR=~/.local/UCSF-Chimera64-1.10.1
    alias chimera="$CHIMERADIR"/bin/chimera
    chimeracli() { chimera --nogui --silent --script "${*}"; }

*Running Chimera from CLI requires the same three parameters all the time: `chimera --nogui --silent --script "<the commands to run>"`. That last function will avoid that repetitive task.*

4 - Save it and apply the changes with `source ~/.bashrc`. This way, you can run Chimera from the console by typing `chimera`. If you want to run a script, simply type:

    chimeracli <script.py> <arg1> <arg2> <...>

**Install Python packages from source**

UCSF Chimera uses its own copy of Python 2.7, so the packages must be installed with that one, and not the system's. For each package, follow these instructions:

1 - Download the source of the package, usually as a *.tar.gz file, and extract it to a new directory with `tar -zxvf <thepackage>.tar.gz`. Once completed, `cd` into that directory.

2 - Instead of using the regular `python setup.py install`, you will install it with Chimera's Python (and our wonderful alias!): 

    chimeracli setup.py install

3 - The bundled Python will copy the contents to the needed folders and precompile the files for you. And that's it!

**Install using pip**

If you are going to be installing a lot of packages, maybe it's a good idea to set up a copy of `pip`, which will resolve dependencies and so on.

1 - `pip` must be installed with `easy_install`, which needs to be installed from source. Download the tar.gz from [here](https://pypi.python.org/pypi/setuptools) (scroll to the bottom) and follow the previous instructions in *Install Python Packages*.

2 - Install `pip` with `easy_install`: 

    chimeracli "$CHIMERADIR"/bin/easy_install pip

3 - Finally, you can run pip installations with `chimeracli "$CHIMERADIR"/bin/pip install <your_package>`, but it's so cumbersome it deserves another alias. In the same fashion, open `~/.bashrc` and add this line at the end: 

    chimerapip() { chimeracli "$CHIMERADIR"/bin/pip "${*}"; }

4 - `pip` supports list of packages, so to install all the dependencies of GAUDI at once, just run these command: 

    chimerapip install deap pyyaml repoze.lru

## Running a GAUDI job
You only have to run `base.py <inputfile>.gaudi` with Chimera's Python. You have already done this to install the packages!

    cd /path/to/gaudi/
    chimeracli base.py /path/to/input/file.gaudi