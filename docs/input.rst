Input files
===========

GAUDI uses YAML-formatted files for both input and output files. YAML is a human-readable serialization format, already implemented in a broad range of languages. Input files use ``gaudi-input`` as file extension, and must contain these five sections:

- **default**. Project options. Configure it to your liking
- **ga**. Genetic algorithm configuration. Normally, you don't have to touch this.
- **similarity**. The similarity function to compare potentially redundant solutions.
- **genes**. List of descriptors used to define an individual
- **objectives**. The list of functions that will evaluate your individuals.
  
You can check some sample input files in the ``examples`` directory.

How to create GAUDI input files
-------------------------------

If you don't mind installing `GAUDInspect <https://bitbucket.org/jrgp/gaudinspect>`_, it provides a full GUI to create GAUDI input files, step by step. Add genes and objectives, configure the paths, number of generations and population size, and run it. Simple and easy.

However, you can also edit them manually, since they are just plain text files. Create a copy of one of the examples and edit them to your convenience. Check the API documentation of :ref:`api.gaudi.genes` and :ref:`api.gaudi.objectives` to find a list of available components.

.. todo::

  * Tutorial: Creating a new input file from scratch