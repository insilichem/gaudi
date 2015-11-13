#!/usr/bin/python2

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
Experimental wrapper to use Chimera's Python directly, instead of a
second subprocess instance inside an already running Python interpreter.
"""

import os
import sys


def chimera_env():
    os.environ['CHIMERA'] = CHIMERA = guess_chimera_path()
    os.environ['CHIMERAPATH'] = CHIMERA
    sys.path.insert(0, os.path.join(CHIMERA, 'share'))
    sys.path.insert(0, os.path.join(CHIMERA, 'bin'))


def chimera_init():
    import chimeraInit
    basepath = os.path.dirname(os.path.abspath(sys.argv[0]))
    chimera_pass = os.path.join(basepath, '__pass__.py')
    chimeraInit.init(['GAUDI', '--script', chimera_pass], nogui=True,
                     eventloop=False, exitonquit=False)


def gaudi_init():
    from gaudi.cli.gaudi_cli import cli
    cli(prog_name='gaudi')


def guess_chimera_path():
    # WINDOWS
    if sys.platform.startswith('win'):
        for programfiles in ('PROGRAMFILES', 'PROGRAMFILES(X86)', 'PROGRAMW6432'):
            paths = [os.path.join(os.environ[programfiles], d)
                     for d in os.listdir(os.environ[programfiles])
                     if d.startswith('Chimera')]
            if paths:
                break
        else:
            paths = [raw_input('Chimera could not be found. Please, type its location manually.\n'
                               'It should be something like C:/Program Files/Chimera 1.10.1\n')]

        paths.sort()
        return paths[-1]  # get latest version available, if more than one exists
    else:
        return subprocess.check_output(['chimera', '--root']).decode('utf-8')


def main():
    chimera_env()
    chimera_init()
    gaudi_init()

if "__main__" == __name__:
    main()
