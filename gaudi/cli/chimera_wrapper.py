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

import os
import sys
import subprocess
import gaudi


def chimera(verbose=False):
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
    chimera(verbose=True)


def guess_windows_chimera():
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
