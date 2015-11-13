#!/usr/bin/python

##############
# GAUDIasm: Genetic Algorithms for Universal
# Design Inference and Atomic Scale Modeling
# Authors:  Jaime Rodriguez-Guerra Pedregal
#            <jaime.rodriguezguerra@uab.cat>
#           Jean-Didier Marechal
#            <jeandidier.marechal@uab.cat>
# Web: https://bitbucket.org/jrgp/gaudi
##############

"""
A wrapper around Chimera script launcher to use a single keyword as binary.

That way, instead of doing chimeracli ``/path/to/gaudi_cli.py``, we can type ``gaudi``
and let setuptools entry_points figure out the rest.

.. note ::
    
    This approach uses subprocess, so theoretically, we could use a queue to launch several
    instances in parallel and simulate multiprocessing. However, there are more adequate 
    strategies we will try first.
"""

import os
import sys
import subprocess
import gaudi


def chimera(verbose=False):
    """
    Find Chimera binary and launch gaudi_cli.py script behind the scenes

    Parameters
    ----------
    verbose : bool, optional
        If True, disable ``--silent`` flag and let Chimera print all it wants
    """
    # Find chimera binary location
    if sys.platform.startswith('linux'):
        chimera = ['chimera']
    elif sys.platform.startswith('win'):
        chimera = [guess_windows_chimera()]
    else:
        sys.exit("ERROR: Platform not supported.")

    # Find GAUDI
    gaudicli = os.path.join(gaudi.__path__[0], 'cli', 'gaudi_cli.py')

    # Prepare command
    script = ' '.join([gaudicli] + sys.argv[1:])
    if verbose:
        command = chimera + ['--nogui', '--debug', '--script', script]
    else:
        command = chimera + ['--nogui', '--debug', '--silent', '--script', script]

    # Launch with clean exit
    sys.exit(subprocess.call(command))


def chimera_verbose():
    """
    Alias to Chimera call with ``--silent`` disabled
    """
    chimera(verbose=True)


def guess_windows_chimera():
    """
    Try to find the Chimera binary in Windows, where the installer does not provide
    a binary in $PATH. It traverses Program Files searching for Chimera installations
    and takes the most recent one.

    .. todo ::

        Untested!
    """
    chimera_paths = [d for d in os.listdir(os.environ['PROGRAMFILES'])
                     if d.startswith('Chimera')]
    if not chimera_paths:  # try with 32 bit
        chimera_paths = [d for d in os.listdir(os.environ['PROGRAMFILES(X86)'])
                         if d.startswith('Chimera')]
    if not chimera_paths:
        user_provided = input('Please, specify the installation folder of UCSF Chimera')
        if os.path.isdir(user_provided):
            chimera_paths = [user_provided]
        else:
            sys.exit("ERROR: Provided directory does not exist.")

    chimera_paths.sort()
    latest = chimera_paths[-1]  # get latest version available, if more than one exists
    return os.path.join(latest, 'bin', 'chimera.exe')
