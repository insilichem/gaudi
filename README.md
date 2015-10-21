# GAUDI_asm Project
GAUDIasm: Genetic Algorithms for Universal Design Inference and Atomic Scale Modeling provides a novel method in design essays by combining several chemostructural criteria, along with energy-like optimization.

## Features
**True multi-objective optimization**

* Feel free to optimize H bonds, hydrophobic interactions, desolvation effects, distances between given sets of atoms, rotamers and more.

**Force-field-less approach**

* Metal complexes are more than welcome.

**Unprecedented customizability**

* Each objective can be switched off upon demand, and multi-instantiated.


## Install Chimera and set aliases in Linux

1 - Download the [latest stable copy of UCSF Chimera](http://www.cgl.ucsf.edu/chimera/download.html) and install it with:

    chmod +x chimera-*.bin && ./chimera-*.bin

2 - Locate the installation directory. The default location is `~/.local/UCSF-Chimera<arch>-<version>`, so adapt it for your version. For example, the binary of a 64-bit Chimera v1.10.1 will be at `~/.local/UCSF-Chimera64-1.10.1`. Write it down for the next step.

3 - Now, create some useful bash aliases. Open `~/.bashrc` with your favourite editor and add these lines at the end of the file. Edit `CHIMERADIR` to suit your installation, but *do not* edit the remaining lines.

    CHIMERADIR=~/.local/UCSF-Chimera<arch>-<version>
    alias chimera="$CHIMERADIR"/bin/chimera
    chimeracli() { chimera --nogui --silent --script "${*}"; }
    chimerapip() { chimeracli "$CHIMERADIR"/bin/pip "${*}"; }


4 - Save it and apply the changes with `source ~/.bashrc`. 

## Install Chimera and set aliases in Windows

1 - Download the [latest stable copy of UCSF Chimera](http://www.cgl.ucsf.edu/chimera/download.html) and open up the .exe to follow the wizard.

2 - Locate the installation directory. The default location is `C:/Program Files/UCSF-Chimera<arch>-<version>`, so adapt it for your version. For example, the binary of a 64-bit Chimera v1.10.1 will be at `C:/Program Files/UCSF-Chimera64-1.10.1`. Write it down for the next step.

3 - Open the Start menu, type `powershell`, and right-click on the first result. Choose `Run as administrator` and check if you already have a profile with the command `Test-Path $Profile`. If the output is `True` go to step 5.

4 - Create a profile with the command:

    New-Item –Path $Profile –Type File –Force

5 - Edit the profile with `notepad $profile` and add these lines, using the proper value for `$CHIMERADIR`.

    $CHIMERADIR = "${env:ProgramFiles}\Chimera 1.10.1"
    Set-Alias chimera "$CHIMERADIR\bin\chimera.exe"
    Set-Alias chimerapip "$CHIMERADIR\bin\Scripts\pip.exe"
    Set-Alias chimerapy "$CHIMERADIR\bin\python.exe"
    function chimeracli { chimera --nogui --debug --script  "$args" }

6 - To apply these changes, run `. $profile`.

## Install GAUDI (platform independent)

With the aliases set above, you can run Chimera from the console by typing `chimera`. Also, if you want to run a Python script, simply type:

    chimeracli <script.py> <arg1> <arg2> <...>

We will use these to set up GAUDI.

1 - You may have noticed we have included an alias called `chimerapip`, which will handle the installation of new Python packages *inside* Chimera. However, Chimera does not include `pip` by default, so you must install that prior to any other package. In order to do so, download `get-pip.py` from [here](https://bootstrap.pypa.io/get-pip.py) and run it with:

    chimeracli /path/to/get-pip.py

2 - Finally, you can run pip installations with `chimerapip`. Download the latest stable GAUDI wheel from [here](https://bitbucket.org/jrgp/gaudi/downloads) and type:

    chimerapip install /path/to/downloaded/gaudi-version.whl


## Running a GAUDI job

You only have to run `launch.py <inputfile>.in.gaudi` with Chimera's Python. Ie:

    chimeracli /path/to/gaudi/scripts/launch.py /path/to/input/file.in.gaudi
