.. GaudiMM: Genetic Algorithms with Unrestricted
   Descriptors for Intuitive Molecular Modeling
   
   http://bitbucket.org/insilichem/gaudi
  
   Copyright 2017 Jaime Rodriguez-Guerra, Jean-Didier Marechal
   
   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at
   
        http://www.apache.org/licenses/LICENSE-2.0
   
   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

GaudiMM
=======

GaudiMM: Genetic Algorithms with Unrestricted Descriptors for Intuitive Molecular Modeling provides a novel method in design essays by combining several chemostructural criteria, along with energy-like optimization.

Features
--------

**True multi-objective optimization**

- Feel free to optimize H bonds, hydrophobic interactions, desolvation effects, distances between given sets of atoms, rotamers and more.

**Force-field-less approach**

- Metal complexes are more than welcome.

**Unprecedented customizability**

- Each objective can be switched off upon demand, and multi-instantiated.


Documentation
-------------

Quick installation:

1 - Download the `latest stable copy of UCSF Chimera <http://www.cgl.ucsf.edu/chimera/download.html>`_ and install it with:

::

  chmod +x chimera-*.bin && sudo ./chimera-*.bin

2 - Install `Miniconda Python 2.7 Distribution <http://conda.pydata.org/miniconda.html>`_ for your platform and install it with:

::

  bash Miniconda2*.sh

3 - Grab the `environment.yml <https://bitbucket.org/insilichem/gaudi/raw/HEAD/environment.yml>`_ file and create the GAUDI environment with:

::

  conda env create -f /path/to/downloaded/environment.yml

4 - Activate the new environment as proposed:

::

  source activate gaudi

5 - Run it!

::

  gaudi

You can also check the (outdated) docs `here <docs/>`_.

Bug reports and contact
-----------------------

Please, use the `issues page <https://bitbucket.org/jrgp/gaudi/issues>`_ of our `Bitbucket repo <https://bitbucket.org/jrgp/gaudi>`_. You can drop me a message at `jaime@klingon.uab.cat <mailto:jaime@klingon.uab.cat>`_ too.