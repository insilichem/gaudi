How to install
==============

GAUDI is written in Python 2.7 and has a main dependency: `UCSF Chimera <https://www.cgl.ucsf.edu/chimera/>`_, *"a highly extensible program for interactive visualization and analysis of molecular structures and related data, including density maps, supramolecular assemblies, sequence alignments, docking results, trajectories, and conformational ensembles"*.

UCSF Chimera has its own embedded Python interpreter, which is heavily modified to achieve some custom behaviour. Using Chimera without that custom Python is hard and troublesome, so we better stick to it. Unfortunately, Chimera does not use any standard environment isolation... which means we have to resort to either Bash aliases, or either dark path magic to make it work easily.

The following steps contain the needed information to install both Chimera and GAUDI in your computer.
  
Install Chimera and set aliases in Linux
----------------------------------------

1 - Download the `latest stable copy of UCSF Chimera <http://www.cgl.ucsf.edu/chimera/download.html>`_ and install it with:

.. code-block:: console

    chmod +x chimera-*.bin && sudo ./chimera-*.bin

2 - Create the desktop and menu shortcuts if you want so. When prompted about creating a symbolic link, please **do it**. Normally, ``/usr/bin`` is fine. (That's why you needed that ``sudo``). If you don't have ``sudo`` permissions, ask your IT admin to create that symlink for you.

3 - Now, create some useful bash aliases. Open ``~/.bashrc`` with your favourite editor and paste these lines at the end of the file.

.. code-block:: bash

    chimeracli() { chimera --nogui --silent --script "${*}"; }
    chimerapip() { chimeracli `chimera --root`/bin/pip "${*}"; }
    chimerapy()  { `chimera --root`/bin/python2.7 ${*}; }


4 - Save it and apply the changes with ``source ~/.bashrc``. 

Install Chimera and set aliases in Windows
------------------------------------------

1 - Download the `latest stable copy of UCSF Chimera <http://www.cgl.ucsf.edu/chimera/download.html>`_ and open up the .exe to follow the wizard.

2 - Locate the installation directory. The default location is ``C:/Program Files/UCSF-Chimera<arch>-<version>``, so adapt it for your version. For example, the binary of a 64-bit Chimera v1.10.1 will be at ``C:/Program Files/UCSF-Chimera64-1.10.1``. Write it down for the next step.

3 - Open the Start menu, type ``powershell``, and right-click on the first result. Choose ``Run as administrator`` and check if you already have a profile with the command ``Test-Path $Profile``. If the output is ``True`` go to step 5.

4 - Create a profile with the command:

.. code-block:: console

    New-Item –Path $Profile –Type File –Force

5 - Edit the profile with ``notepad $profile`` and add these lines, using the proper value for ``$CHIMERADIR``.

.. code-block:: console

    $CHIMERADIR = "${env:ProgramFiles}\Chimera 1.10.1"
    Set-Alias chimera "$CHIMERADIR\bin\chimera.exe"
    Set-Alias chimerapip "$CHIMERADIR\bin\Scripts\pip.exe"
    Set-Alias chimerapy "$CHIMERADIR\bin\python.exe"
    function chimeracli { chimera --nogui --debug --script  "$args" }

6 - To apply these changes, run ``. $profile``.

Install GAUDI (platform independent)
------------------------------------

With the aliases set above, you can run Chimera from the console by typing ``chimera``. Also, if you want to run a Python script, simply type:

.. code-block:: console

    chimeracli <script.py> <arg1> <arg2> <...>

We will use these to set up GAUDI.

1 - You may have noticed we have included an alias called ``chimerapip``, which will handle the installation of new Python packages *inside* Chimera. However, Chimera does not include ``pip`` (the installer utility) by default, so you must install that prior to any other package. In order to do so, download ``get-pip.py`` from `here <https://bootstrap.pypa.io/get-pip.py>`_ and run it with:

.. code-block:: console

    chimeracli /path/to/downloaded/get-pip.py

2 - Finally, you can run pip installations with ``chimerapip``. For GAUDI, use this (long) command:

.. code-block:: console

    chimerapip install --extra-index-url http://klingon.uab.cat/repo/jaime/pip --trusted-host klingon.uab.cat --allow-unverified gaudi gaudi



Check everything is OK
----------------------

If everything went OK, you will have a ``gaudi`` binary along Chimera binaries. Link it to somewhere in your ``$PATH``. In Linux, it's something like:

.. code-block:: console
    
    # Linux
    sudo ln -s `chimera --root`/bin/gaudi /usr/local/bin/gaudi

For Windows, you have to open an administrator ``cmd.exe``. You will find a shortcut for that if you press ``Win+X``. Then, type:

.. code-block:: console

    # Windows Vista/7+
    mklink $CHIMERADIR/bin/gaudi.exe C:/WINDOWS/gaudi.exe

Now, if you type ``gaudi``, you will get the usage screen:

.. code-block:: console

    Usage: gaudi [OPTIONS] COMMAND [ARGS]...

      GAUDI: Genetic Algorithms for Universal Design Inference

      By Jaime Rodríguez-Guerra and Jean-Didier Maréchal.
      https://bitbucket.org/jrgp/gaudi

    Options:
      --version   Show the version and exit.
      -h, --help  Show this message and exit.

    Commands:
      benchmark  Performs the same essay over a dataset.
      prepare    Create or edit a GAUDI input file.
      rmsd       Calculate RMSD of results vs reference.
      run        Launch a GAUDI input file.
      view       Analyze the results in a GAUDI output file.


However, if that doesn't work, there is a manual method you can alias in your ``.bashrc``:

.. code-block:: console
    
    gaudi() { chimeracli `chimera --root`/lib/python2.7/site-packages/gaudi/cli/gaudi_cli.py ${*}; }


... or PowerShell profile:

.. code-block:: console
    
    function gaudi { chimeracli $CHIMERADIR/lib/python2.7/site-packages/gaudi/cli/gaudi_cli.py $args }