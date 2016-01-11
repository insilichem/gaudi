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

from __future__ import print_function
import os
import sys
import subprocess


def chimera_env(force_ipython=False):
    """
    Patch current environment variables so Chimera can start up and we can import its modules
    """
    if 'CHIMERA' not in os.environ:
        CHIMERA = os.environ['CHIMERA'] = os.environ['PYTHONHOME'] = guess_chimera_path()

        # PYTHONPATH defines the locations where Python should look for packages and modules
        # The sys.path corresponds to the system python environment path (conda, probably)
        os.environ['PYTHONPATH'] = ":".join([os.path.join(CHIMERA, 'bin'),
                                             os.path.join(CHIMERA, 'share'),
                                             os.path.join(CHIMERA, 'lib')]
                                            + sys.path)

        # Load Chimera libraries
        CHIMERALIB = os.path.join(CHIMERA, 'lib')
        if sys.platform == 'win32':
            os.environ['PATH'] += ":" + CHIMERALIB
        elif sys.platform == 'darwin':
            OLDLIB = os.environ.get('DYLD_LIBRARY_PATH', '')
            os.environ['DYLD_LIBRARY_PATH'] = ':'.join([CHIMERALIB, OLDLIB])
        else:
            OLDLIB = os.environ.get('LD_LIBRARY_PATH', '')
            os.environ['LD_LIBRARY_PATH'] = ':'.join([CHIMERALIB, OLDLIB])

        # Reload Python with modified environment variables
        if inside_ipython() or force_ipython:
            executable = os.path.join(os.path.dirname(sys.executable), 'ipython')
            os.environ['TERM'] = "xterm-256color"
        else:
            executable = sys.executable

        if interactive_mode() or force_ipython:
            sys.argv.insert(0, '-i')

        try:
            os.execve(executable, [executable] + sys.argv, os.environ)
        except OSError:
            sys.exit("ERROR: {} not found.\n"
                     "Have you installed IPython in your environment? Use:"
                     "    conda install ipython".format(executable))


def chimera_init():
    """
    Bypass script loading and initialize Chimera in nogui mode.
    """
    import chimeraInit
    from chimera import registration
    from Midas import midas_text

    # bypass ReadStdin extension & registration
    old_doRunScript = midas_text.doRunScript
    registration.checkRegistration = midas_text.doRunScript = lambda *args: None
    chimeraInit.init(['GAUDI', '--nostatus', '--silent', '--script', ''], nogui=True,
                     eventloop=False, exitonquit=False)
    midas_text.doRunScript = old_doRunScript
    del registration, chimeraInit, midas_text


def chimera_test():
    """
    Test if Chimera can be imported
    """
    try:
        import chimera
    except ImportError:
        print('Warning: Could not import Chimera...')
    else:
        print()
        print('-'*70)
        msg = 'UCSF Chimera is now enabled in this interpreter'
        padding = ((70-len(msg))/2)-1
        print(' '*padding, msg, ' '*padding)
        print('-'*70)


def gaudi_init():
    """
    Launch GAUDI
    """
    from gaudi.cli.gaudi_cli import cli
    cli(prog_name='gaudi')


def guess_chimera_path():
    """
    Try to guess Chimera installation path
    """
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


def inside_ipython():
    try:
        __IPYTHON__
        return True
    except NameError:  # IPython not in use
        return False


def interactive_mode():
    return (hasattr(sys, 'ps1') and sys.ps1) or sys.flags.interactive


def main():
    chimera_env()
    chimera_init()


def main_ipython():
    chimera_env(force_ipython=True)
    chimera_init()
    chimera_test()


def main_with_gaudi():
    main()
    gaudi_init()


if "__main__" == __name__:
    main()
    chimera_test()
