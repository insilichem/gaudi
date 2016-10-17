#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function
from setuptools import setup, find_packages
import os
import io

import versioneer

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
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    url=gaudi.__url__,
    license='Free For Educational Use',
    author=gaudi.__author__,
    author_email='jaime.rogue@gmail.com',
    description=gaudi.__description__,
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
    # dependencies are handled with conda-recipe/meta.yaml
    # check that file if you used setup.py manually
    # install_requires='deap click repoze.lru bunch voluptuous pyyaml'.split(),
    entry_points='''
        [console_scripts]
        gaudi=gaudi.cli.gaudi_cli:cli
        '''
    ,
)
