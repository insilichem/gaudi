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

.. _primer:

============================
GaudiMM basics (start here!)
============================

The main strength of GaudiMM is, possibly, its flexible approach to solve molecular modeling problems. However, flexibility does not come without a price: the learning curve can be somewhat steep for beginners. In this tutorial we will try to smooth that curve out.

GaudiMM makes a big effort to clearly separate *exploration* from *evaluation*. As such, both stages are configured in different sections of the input file: *genes* and *objectives*.

Exploration and genes
=====================

Most molecular modeling problems can be explained using search space analogies. This is, given a problem, one can identify the dimensions it must explore to find adequate solutions for that problem. That n-dimensional space can be comprised of one or more (bio)chemical structures and their respective spatial variability. Modifying the topology and/or composition of such structures would be equivalent to moving along the (bio)chemical axes of such search space, while modifying their atomic positions would be equivalent to moving along the cartesian and/or internal coordinates axes. Thus, visiting different regions of the whole search space means creating new candidate solutions to the problem... We will call this stage *exploration*, and it is detailed in the genes section of the input file. If you want to know why they are called genes, check the *Genetic Algorithms Primer* in the :ref:`developers`.

The more important gene is, almost always, the :class:`gaudi.genes.molecule.Molecule` gene. This is used to load (bio)chemical structures, such as proteins, peptides and small cofactors. Normally they are straight PDB or Mol2 files, but can also be more complex constructs. As of now, it is the only one capable of providing topologies.
The other genes are responsible of applying changes to the coordinates of that initial structure. Some examples include torsion angles of rotatable bonds (:class:`gaudi.genes.torsion.Torsion` gene), translating and rotating a Molecule in the 3D space (:class:`gaudi.genes.search.Search` gene), or importing the atomic positions of a molecular dynamics trajectory frame (:class:`gaudi.genes.trajectory.Trajectory` gene).

Evaluation and objectives
=========================

However, there is no guarantee that, given a new position on such search space, that hypothetical solution is actually good enough to solve the problem. To assess the quality of a random solution, it must be *evaluated* with some criteria or *objectives*. 
All objectives have something in common: they take a candidate solution (ie, a point of the search space) and measure some property to return a numeric score. Some examples include the distance between two given atoms (:class:`gaudi.objectives.distance.Distance` objective), the steric clashes of the molecules (:class:`gaudi.objectives.contacts.Contacts` objective) or the forcefield energy of the system (:class:`gaudi.objectives.energy.Energy` objective).

Examples of application
=======================

Now we got the basics covered, we can start to think of how we could solve some standard (and not so standard) molecular modeling problems with GaudiMM.

Docking
-------

Docking problems devote to finding the correct orientation and position of a small molecule (the ligand) within the cavity of a bigger one (normally, a protein). For the exploration stage, at least three genes must be defined:

- A :class:`gaudi.genes.molecule.Molecule` gene to load the file corresponding to the protein structure.
- Another :class:`gaudi.genes.molecule.Molecule` gene to load the ligand itself.
- A :class:`gaudi.genes.search.Search` gene to translate and rotate the ligand around the protein surroundings.
  
This is all you need to move a rigid ligand around a protein. However, most of the time you want to implement some kind of internal **flexibility**, such as torsions in the ligand or some residues sidechains. To do that, add these to your ``genes`` section:

- A :class:`gaudi.genes.torsion.Torsion` gene with the ligand as target to implement dihedral angle torsions along the rotatable bonds of the ligand.
- A :class:`gaudi.genes.rotamers.Rotamers` gene with some residues of the protein as target. This will randomly modify the dihedral torsions of the sidechains of the specified residues according to the angles provided by the rotamer library (defaults to Dunbrack's).
  
There is also a :class:`gaudi.genes.mutamers.Mutamers` gene that allows you to mutate residues in addition to rotameric exploration. This will provide some traversal along the biochemical axis of your protein, but careful, because the current implementation is RAM hungry!

The *evaluation* stage can be comprised of several objectives, but normally you'd want to:

- Minimize steric clashes with :class:`gaudi.objectives.contacts.Contacts`.
- Maximize hydrophobic patches with another :class:`gaudi.objectives.contacts.Contacts` entry.
- Use a proper docking scoring function, such as :class:`gaudi.objectives.ligscore.LigScore` or :class:`gaudi.objectives.dsx.DSX`. Depending on the one you choose, it should be minimized (the usual case) or maximized.

Covalent docking
----------------

Not all docking studies consider a free ligand. Sometimes the ligand is anchored to some part of the protein. While there's no specific gene to implement a covalent bond (for now, at least), you can mimick it with a Search gene with ``radius=0`` and ``rotate=True``. The origin of this null search sphere should be set to the atom that is part of the covalent bond, so you would have add an extra atom in one of the molecules. Also, we recommend setting an :class:`gaudi.objectives.angle.Angle` objective between the involved atoms in the covalent bond so that the resulting rotation matches the expected geometry of the new bond (109.5º for *sp3*, 120º for *sp2*, 180º for *sp1*).

Naked metal ions docking
------------------------

If instead a small organic molecule you choose to load a naked metal ion... would it work? Well, in most docking programs probably not, but with GaudiMM the answer is... it depends!

It depends on the objectives you choose. Since we ship a :class:`gaudi.objectives.coordination.Coordination` objective that can deal with coordination geometries, you can dock naked metal ions in any protein or peptide.

Competitive docking
-------------------

Since you can instantiate as many :class:`gaudi.genes.molecule.Molecule` genes as you want, nothing prevents you from adding more ligand molecules at the same time.  They will compete to find its place in the protein(s). Just remember to replicate any dependent genes (:class:`gaudi.genes.search.Search`, :class:`gaudi.genes.torsion.Torsion` and so on) accordingly.
Also, why not two proteins which the ligand should choose from (but it's true there are finer ways to assess that)?

Hacking Molecule genes for complex studies
------------------------------------------

In addition to loading molecular structure files, the :class:`gaudi.genes.molecule.Molecule` does a couple of extra things. The path parameter can be set to two different values:

- The path to a PDBor mol2 file (or any format that Chimera can open). This is the standard behaviour. It will load the structure and that's it.

- The path to a directory, whose contents determine the final behaviour:

    A) If the directory contains molecule files, GaudiMM will choose one of them randomly for each case.

    B) If the directory contains subdirectories which, in turn, contain molecules files, GaudiMM will sort those subdirectories by name and then pick one molecule from each, in that order. The chosen molecules will be chained linearly as specified in the accompanying ``*.attr`` files. Perfect for drug design, since you can fill a couple of directories with different alkanes to link a cofactor

Case **A** allows GaudiMM to deal test a library of compounds against certain criteria: ie, virtual screening!

Case **B** makes drug design studies possible. If you want to test different linker molecules to anchor a cofactor to a protein, you can use a two-directories scheme: one directory would contain the aforementioned linkers, and the other one just the cofactor (it's fine it there's only one molecule in one directory).

Peptide folding & Conformational analysis
=========================================

Loading an unfolded peptide with :class:`gaudi.genes.molecule.Molecule` is pretty easy. Just specify the path. Then, if you apply a :class:`gaudi.genes.torsion.Torsion` gene on such peptide, but considering only the bonds involving an alpha-carbon, you can effectively explore the conformational space of its folding. Use :class:`gaudi.genes.rotamers.Rotamers` to implement sidechain flexibility.

For the evaluation, you can use the :class:`gaudi.objectives.energy.Energy` objective to assess the forcefield energy as provided by OpenMM GPU calculations. This type of study could be considered some sort of highly explosive metadynamics.

Homology modeling
-----------------

In the same fashion, given an unfolded peptided or protein segment, you could apply the same Torsion scheme to explore different conformations. A hypothetical objective could be devised to calculate the RMS deviation of the new fragment and a reference one, which then the GA would minimize to obtain a similar structure. 

Conformational analysis
-----------------------

Another possible variation of this scheme is to impose geometric guides as additional objectives, like :class:`gaudi.objectives.distance.Distance`, :class:`gaudi.objectives.angle.Angle` or :class:`gaudi.objectives.coordination.Coordination`. The resulting calculation could be regarded as a restrained conformational analysis, very useful for finding initial structures of unparametrized small molecules you want to study with higher levels of theory, such as QM.

Trajectory analysis
-------------------

GaudiMM features a :class:`gaudi.genes.trajectory.Trajectory` gene capable of importing frames of MD movies and apply them to the corresponding :class:`gaudi.genes.molecule.Molecule` instance. This way, all the objectives are available as trajectory analysis tools, maybe not present in other specific software. For example, finding coordination geometries of a given metal ion along a molecular mechanics simulation.