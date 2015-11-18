#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function
from setuptools import setup, find_packages
import os
import io

import gaudi

here = os.path.abspath(os.path.dirname(__file__))


def read(*filenames, **kwargs):
    encoding = kwargs.get('encoding', 'utf-8')
    sep = kwargs.get('sep', '\n')
    buf = []
    for filename in filenames:
        with io.open(os.path.join(here, filename), encoding=encoding) as f:
            buf.append(f.read())
    return sep.join(buf)

long_description = read('README.rst')

setup(
    name='gaudi',
    version=gaudi.__version__,
    url='https://bitbucket.org/jrgp/gaudi',
    license='Free For Educational Use',
    author='Jaime Rodriguez-Guerra Pedregal',
    author_email='jaime.rogue@gmail.com',
    description='A full GUI for launching GAUDI jobs, '
                'analyzing their progress, and examining their results.',
    long_description=long_description,
    packages=find_packages(),
    include_package_data=True,
    platforms='any',
    classifiers=[
        'Programming Language :: Python :: 2.7',
        'Development Status :: 2 - Pre-Alpha',
        'Natural Language :: English',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: Free For Educational Use',
        'Operating System :: OS Independent',
        'Topic :: Scientific/Engineering :: Chemistry',
    ],
    install_requires=[
        'PyYaml',
        'deap',
        'repoze.lru',
        'click'
    ],
    entry_points='''
        [console_scripts]
        gaudi=gaudi.cli.chimera_wrapper:chimera
        gaudiv=gaudi.cli.chimera_wrapper:chimera_verbose
    ''',
    # gaudi=gaudi.cli.gaudi_cli:cli
)
