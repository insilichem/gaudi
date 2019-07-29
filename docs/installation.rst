.. GaudiMM: Genetic Algorithms with Unrestricted
   Descriptors for Intuitive Molecular Modeling

   https://github.com/insilichem/gaudi

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

==============
How to install
==============

Quick steps:

1 - Download the `latest stable copy of UCSF Chimera <http://www.cgl.ucsf.edu/chimera/download.html>`_ and install it with:

::

  chmod +x chimera-*.bin && sudo ./chimera-*.bin

2 - Install `Miniconda Python 2.7 Distribution <http://conda.pydata.org/miniconda.html>`_ for your platform and install it with:

::

  bash Miniconda2*.sh

3 - Install ``gaudi`` with ``conda`` in a new environment called ``insilichem`` (or whatever name you prefer after the ``-n`` flag), using these custom channels (``-c`` flags):

::

  conda create -n insilichem -c omnia -c salilab -c insilichem -c anaconda -c bioconda gaudi


4 - Activate the new environment as proposed:

::

  source activate insilichem

5 - Run it!

::

  gaudi


Check everything is OK
======================

If everything went OK, you will get the usage screen:

.. code-block:: console

    Usage: gaudi [OPTIONS] COMMAND [ARGS]...

      GaudiMM: Genetic Algorithms with Unrestricted Descriptors for Intuitive
      Molecular Modeling

      (C) 2017, InsiliChem
      https://github.com/insilichem/gaudi

    Options:
      --version   Show the version and exit.
      -h, --help  Show this message and exit.

    Commands:
      prepare  Create or edit a GAUDI input file.
      run      Launch a GAUDI input file.
      view     Analyze the results in a GAUDI output file.
