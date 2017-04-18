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


.. _tutorial:

========================================
Tutorial: Your first GaudiMM calculation
========================================

Running a GaudiMM calculation is easy. All you need is a simple command:

.. code-block:: console

    gaudi run input_file.yaml

The question is... how do I create that ``input_file.yaml``.

How to create an input file from scratch
----------------------------------------

Input files in GaudiMM are formatted with YAML, which follows readable conventions that are still parseable by computers. To easily edit YAML files, we recommend using a text editor that supports syntax highlighting, such as `Sublime Text 3 <https://www.sublimetext.com/3>`_, `Atom <https://atom.io/>`_ or `Visual Studio Code <https://code.visualstudio.com/>`_. If you use ``.yaml`` as the extension of the input filename, the editor will colorize the text automatically. Else, you can always configure it manually. Check the docs of your editor to do so. If you don't want syntax highlighting it's fine, GaudiMM will work the same.

Input files must contain the sections mentioned above: output, ga, similarity, genes, and objectives. Each of this sections it's defined with a colon and a newline, and its contents will be inside an indented block. For readability, I usually insert a blank line between sections.

.. code-block:: yaml

    output:
        contents of output

    ga:
        contents of ga

    similarity:
        contents of similarity

    genes:
        contents of genes

    objectives:
        contents of objectives

.. note::
    
    The API documentation of :class:`gaudi.parse.Settings` contains the full list of parameters for ``output`` and ``ga`` sections. Programmatically defined default values are always defined in :attr:`gaudi.parse.Settings.default_values`. The appropriate types and whether they are required or not are defined in :attr:`gaudi.parse.Settings.schema`. GaudiMM will check if the submitted values conform to these rules and report any possible mistakes.


The output section
..................

This section governs how the results and reports will be created. Everything is optional, since each key has a default value, but it is preferrable to at least specify the ``name`` of the job (otherwise, it will be set to five random characters), and the ``path`` where the result files will be written. Like this:

.. code-block:: yaml

    output:
        name: some_example
        path: results

.. note::

    All the relative paths in GaudiMM are relative to the location of the input file, not the working directory. To ease this difference, we recommend running the jobs from the same folder where the input file is located. Of course, you can always use absolute paths.

The ga section
..............

This section hosts the parameters of the Genetic Algorithm GaudiMM uses. Unless you know what you are doing, the only values you should modify are ``population`` and ``generations``. To see the appropriate values, refer to :ref:`faq`. For example:

.. code-block:: yaml

    ga:
        population: 200
        generations: 100

The similarity section
......................

This section contains the parameters to the similarity operator, which, given two individuals with the same fitness, whether they can be considered the same solution or not. This section is deliberatedly loose: you define the Python function to call, together with its positional and keyword arguments. 

For the time being, the only similarity function we ship is based on the RMSD of the two structures: :func:`gaudi.similarity.rmsd`. The arguments are which ``Molecule`` genes should be compared and the RMSD threshold to consider whether they are equivalent or not.

.. code-block:: yaml

    similarity:
        module: gaudi.similarity.rmsd
        args: [[Ligand], 1.0]
        kwargs: {}

The genes section
.................

This section describes the components of the exploration stage of the algorithm; ie, the features of each Individual in the population. While the previous sections were dictionaries (this is, a collection key-value pairs), the ``genes``  and ``objectives`` section is actually a list of dictionaries. As a result, you need to specify them like this:

.. code-block:: yaml

    genes: 
        -   name: Protein
            module: gaudi.genes.molecule
            path: /path/to/protein.mol2

        -   name: Torsion
            module: gaudi.genes.torsion
            target: Ligand
            flexibility: 360

Notice the dash ``-`` next to ``name``. This, and the extra indentation, define a list. Each element of this list is a new gene. Each gene must include two compulsory values:

- ``name``. A unique identifier for this gene. If you add two genes with the same name, GaudiMM will complain.
- ``module``. The Python import path to the module that contains the gene. All GaudiMM builtin genes are located at :mod:`gaudi.genes`.
  
All other parameters are determined by the chosen gene. Check the corresponding documentation for each one!

.. note::
    
    | **How do I know which genes to use?**
    | Unless you code a gene of your own to replace it, you will always need one or more :class:`gaudi.genes.molecule.Molecule` genes. Then, choose the flexibility models you want to implement on top of such molecule. Several examples are provided in :ref:`primer`.

The objectives section
......................

Like the genes section, the objectives section is also a list of dictionaries, so they follow the same syntax:

.. code-block:: yaml

    -   name: Clashes
        module: gaudi.objectives.contacts
        which: clashes
        weight: -1.0
        probes: [Ligand]
        radius: 5.0

    -   name: LigScore
        module: gaudi.objectives.ligscore
        weight: -1.0
        proteins: [Protein]
        ligands: [Ligand]
        method: pose

In addition to the required ``name`` and ``module`` parameters, each objective needs a ``weight`` parameter. If set to ``1.0``, the algorithm will maximize the score returned by the objective; if set to ``-1.0``, it will be minimized. Theoretically, any other positive or negative float will work, but stick to the convention of using ``1.0`` or ``-1.0``. 

Any other parameters present in an objective are responsibility of that objective, and are specified in its corresponding documentation.

That's it! Now save it with a memorable filename and run it!

How to run your input file
--------------------------

Let's get back to the beginning of the tutorial: all you need to do is typing:

.. code-block:: console

    gaudi run input_file.yaml

If everything is fine, you'll see the following output in the console:

.. code-block:: console
    
    $> gaudi run input_file.yaml

      .g8"""bgd       db   `7MMF'   `7MF'`7MM"""Yb. `7MMF'                                    
    .dP'     `M      ;MM:    MM       M    MM    `Yb. MM                                      
    dM'       `     ,V^MM.   MM       M    MM     `Mb MM  `7MMpMMMb.pMMMb.  `7MMpMMMb.pMMMb.  
    MM             ,M  `MM   MM       M    MM      MM MM    MM    MM    MM    MM    MM    MM  
    MM.    `7MMF'  AbmmmqMA  MM       M    MM     ,MP MM    MM    MM    MM    MM    MM    MM  
    `Mb.     MM   A'     VML YM.     ,M    MM    ,dP' MM    MM    MM    MM    MM    MM    MM  
      `"bmmmdPY .AMA.   .AMMA.`bmmmmd"'  .JMMmmmdP' .JMML..JMML  JMML  JMML..JMML  JMML  JMML.
    ------------------------------------------------------------------------------------------
    GaudiMM: Genetic Algorithms with Unrestricted Descriptors for Intuitive Molecular Modeling
    2017, InsiliChem · v0.0.2+251.g122cdf0.dirty

    Loaded input input_file.yaml
    Launching job with...
      Genes: Protein, Ligand, Rotamers, Torsion, Search
      Objectives: Clashes, Contacts, HBonds, LigScore

After the first iteration is complete, the realtime report data will kick in:

.. code-block:: console

    gen progress    nevals  speed       eta     avg                                                 std                                                 min                                     max                                              
    0   4.76%       20      1.25 ev/s   0:16:34 [  3.080e+03  -2.690e+02   7.500e-01   1.179e+03]   [  1.027e+03   7.200e+01   8.292e-01   4.964e+02]   [ 584.881 -398.517    0.     144.78 ]   [  4.753e+03  -1.053e+02   3.000e+00   2.066e+03]
    1   9.52%       60      1.28 ev/s   0:16:32 [  2.787e+03  -2.659e+02   1.400e+00   9.142e+02]   [ 1198.699    98.675     1.393   484.309]           [ 415.912 -398.517    0.      16.45 ]   [ 4538.766   -71.562     5.     1854.63 ]        

The first time you see this it might result too confusing, especially if the terminal wraps long lines. Let's describe each tab-separated column:

- **gen**. The current generation.
- **progress**. Percentage of completion of the job. This is estimated with the expected number of operations: :math:`(generations + 1) * lambda\_ * (cxpb + mutpb)`
- **nevals**. Number of evaluations performed in current generation. 
- **speed**. Estimated number of evaluations per second. This does not take into account the time spent in the variation stage.
- **eta**. Estimated time left.
- **avg**. Average of all the fitness values reported by each objective in the current generation. They are listed in the order given in the input file, and also reflected above, after the *Launching job with...* line.
- **std**. Same as avg, but for the standard deviation.
- **max**. The maximum fitness value reported by each objective in the current generation.
- **min**. Same as above, but for the minimum value.

If the setting ``check_every`` in the ``output`` section is greater than zero, GaudiMM will dump the current population every ``check_every`` generations. That way, you can assess the progress visually along the simulation. 

Also, if you feel that the algorithm has progressed enough to satisfy your needs, you can cancel it prematurely with ``Ctrl+C``. GaudiMM will detect the interruption and offer to dump the current state of the simulation: 

.. code-block:: console

    ^C[!] 

    Interruption detected. Write results so far? (y/N):

Answer ``y`` and wait a couple of seconds while GaudiMM writes the results. To analyze them, check the following tutorial:

- :ref:`tutorial-visualization`
