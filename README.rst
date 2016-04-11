GAUDI Project
=============
GAUDIasm: Genetic Algorithms for Universal Design Inference and Atomic Scale Modeling provides a novel method in design essays by combining several chemostructural criteria, along with energy-like optimization.

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