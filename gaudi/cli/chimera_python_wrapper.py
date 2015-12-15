#!/usr/bin/env python2
# -*- coding: utf-8 -*-

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
import subprocess


def chimera_env():
    if 'CHIMERA' not in os.environ:
        CONDA = guess_conda_env_path()
        CHIMERA = os.environ['CHIMERA'] = os.environ['PYTHONHOME'] = guess_chimera_path()

        # PYTHONPATH defines the locations where Python should look for packages and modules
        os.environ['PYTHONPATH'] = ":".join([
            os.path.join(CHIMERA, 'bin'),
            os.path.join(CHIMERA, 'lib'),
            os.path.join(CHIMERA, 'share'),
            os.path.join(CHIMERA, 'lib', 'python2.7', 'site-packages'),
            os.path.join(CONDA, 'lib', 'python2.7', 'site-packages')])

        if sys.platform == 'win32':
            os.environ['PATH'] += ":" + os.path.join(CHIMERA, 'lib')
        elif sys.platform == 'darwin':
            os.environ['DYLD_LIBRARY_PATH'] = '{}:{}'.format(os.environ.get('DYLD_LIBRARY_PATH', ''),
                                                             os.path.join(CHIMERA, 'lib'))
        else:
            os.environ['LD_LIBRARY_PATH'] = '{}:{}'.format(os.environ.get('LD_LIBRARY_PATH', ''),
                                                           os.path.join(CHIMERA, 'lib'))
        # Reload Python with modified environment variables
        os.execve(sys.executable, [sys.executable] + sys.argv, os.environ)


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
        try:
            return subprocess.check_output(['chimera.exe', '--root']).decode('utf-8').strip()
        except (OSError, subprocess.CalledProcessError):
            for programfiles in ('PROGRAMFILES', 'PROGRAMFILES(X86)', 'PROGRAMW6432'):
                paths = [os.path.join(os.environ[programfiles], d)
                         for d in os.listdir(os.environ[programfiles])
                         if os.path.isdir(d) and d.startswith('Chimera')]
                if paths:
                    break
            else:
                paths = [raw_input('Chimera could not be found. Please, type its location manually.\n'
                                   'It should be something like C:/Program Files/Chimera 1.10.1\n')]

            paths.sort()
            return paths[-1]  # get latest version available, if more than one exists
    else:
        try:
            return subprocess.check_output(['chimera', '--root']).decode('utf-8').strip()
        except (OSError, subprocess.CalledProcessError):
            local = os.path.expanduser('~/.local')
            paths = [os.path.join(local, d) for d in os.listdir(local)
                     if os.path.isdir(d) and d.startswith('UCSF-Chimera')]
            if paths:
                paths.sort()
                return paths[-1]
            else:
                path = raw_input('Chimera could not be found. Please, type its location manually.\n'
                                 'It should be something like /home/user/.local/UCSF-Chimera-1.10.1\n')
                return path.strip()


def guess_conda_env_path():
    CONDA = subprocess.check_output(['conda', 'info', '--root']).decode('utf-8').strip()
    gaudi_default = os.path.join(CONDA, 'envs', 'gaudi')
    if os.path.exists(gaudi_default):
        return gaudi_default
    return ''


def main():
    chimera_env()
    chimera_init()
    gaudi_init()

if "__main__" == __name__:
    main()
