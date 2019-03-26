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

.. _faq:

==================
FAQ & Known issues
==================

Known issues
============

Most problems we have faced have to do with the installation and Chimera-Conda incompatibilities. Those are covered in the `pychimera documentation <https://pychimera.readthedocs.io/en/latest/faq.html>`_, so please check that list before raising an issue in this repository!

How can I cite GaudiMM? What licenses apply?
============================================

GaudiMM is scientific software, funded by public research grants. If you make use of GaudiMM in scientific publications, please cite it. It will help measure the impact of our research and future funding!

.. code-block:: latex
   
   @preamble{ " \newcommand{\noop}[1]{} " }
   @article{GaudiMM2017,
      title = {GaudiMM: A Modular Multi-Objective Platform for Molecular Modeling},
      author = {Rodr{\'i}guez-Guerra Pedregal, Jaime and Sciortino, Giuseppe and Guasp, Jordi and Municoy, Mart{\'i} and Mar{\'e}chal, Jean-Didier},
      year = {\noop{2017}submitted},
   }

GaudiMM itself is licensed under `Apache License 2.0 <https://www.apache.org/licenses/LICENSE-2.0.html>`_, but includes work from other developers, whose licenses apply. Please check the ``LICENSE`` file in the root directory for further details.


How many generations / Which population size should I pick?
===========================================================
This depends on the complexity of your search space (related to the number and type of genes in use), and the evaluation power of the chosen objectives. A simple job would work OK with 50 generations and populations of 100 individuals, while more complex ones would require 500 generations for populations of 1000 indivuals. This also depends on the values of mu and lambda parameters.


The output produces a lot of similar results and I want to focus on diversity
=============================================================================
You should play with the cutoff value of the RMSD ``similarity`` section. If the RMSD of two potentially similar solutions is under the cutoff, one of them is discarded. However, this is only applied if the score of two solutions are the same; ie, if there is a draw.

To force draws, one can reduce the number of decimal positions returned by every objective, which is controlled by the ``precision`` parameter. It is setup globally in the ``output`` section, but can also be overriden by any objective. By default, ``precision`` is 3. Also, if you are using the :class:`gaudi.genes.search.Search` gene, make sure its ``precision`` parameter is small enough so you don't lose exploration efforts in too similar positions.


The output produces very different results and I am interested in clustered solutions
=====================================================================================
In this case, one should increase the ``precision`` value (both globally and in the :class:`gaudi.genes.search.Search` gene, if it applies) and reduce the RMSD ``similarity`` cutoff. 


What is the difference between ``atom_names`` and ``atom_types`` in some objectives or genes?
=============================================================================================
``atom_names`` refers to the ``name`` attribute of ``chimera.Atom`` objects, while ``atom_types`` is applied with ``idatmType`` attributes. Normally, the ``name`` attribute is picked directly from the molecule file (PDB, mol2), while the ``idatmType`` is assigned algorithmically by UCSF Chimera.


.. tip::
   
   Any further questions? Feel free to submit your inquiries to our `issues page <https://github.com/insilichem/gaudi/issues>`_!
