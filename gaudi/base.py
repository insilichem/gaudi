#!/usr/bin/python

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
:mod:`gaudi.base` contains the core classes we use to build individuals
(potential solutions of the optimization process).
"""

# Python
from collections import OrderedDict
from zipfile import ZipFile, ZIP_DEFLATED, ZIP_STORED
import os
import pprint
import math
import logging
# Chimera
import chimera
# External dependencies
import deap.base
import yaml
# GAUDI
import gaudi.plugin

pp = pprint.PrettyPrinter(4)

logger = logging.getLogger(__name__)


class Individual(object):

    """
    Base class for `individual` objects that are evaluated by DEAP.

    Each individual is a potential solution. It contains all that is needed
    for an evaluation. With multiprocessing in mind, individuals should be
    self-contained so it can be passed between threads.

    The defined methods are only wrapper calls to the respective methods of each
    gene.

    :cfg:   The full parsed object from the configuration YAML file.
    :cache: A mutable object that can be used to store values across instances.

    .. todo::

        :meth:`write` should use Pickle and just save the whole object, but
        Chimera's inmutable objects (Atoms, Residues, etc) get in the way. A
        workaround may be found if we take a look a the session saving code.
    """

    _CACHE = {}
    _CACHE_OBJ = {}

    def __init__(self, cfg=None, cache=None, **kwargs):
        logger.debug("Creating new individual with id %s", id(self))
        self.genes = OrderedDict()
        self.cfg = cfg
        gaudi.plugin.load_plugins(self.cfg.genes, container=self.genes,
                                  cache=self._CACHE, parent=self,
                                  cxeta=self.cfg.ga.cx_eta,
                                  mteta=self.cfg.ga.cx_eta,
                                  indpb=self.cfg.ga.mut_indpb)
        for g in self.genes.values():
            g.__ready__()

        self.fitness = gaudi.base.Fitness(parent=self, cache=self._CACHE_OBJ)

    def evaluate(self):
        logger.debug("Evaluating individual #%s", id(self))
        self.express()
        score = self.fitness.evaluate()
        self.unexpress()
        return score

    def express(self):
        """
        Express genes in this environment. Very much like 'compiling' the
        individual to a chimera.Molecule.
        """
        for name, gene in self.genes.items():
            logger.debug("Expressing gene %s with allele\n%s",
                         name, pp.pformat(gene.allele))
            gene.express()

    def unexpress(self):
        """
        Undo .express()
        """
        for gene in reversed(self.genes.values()):
            logger.debug("Reverting expression of gene %s", gene.name)
            gene.unexpress()

    def mate(self, other):
        for gene in self.genes.values():
            gene.mate(other.genes[gene.name])
        logger.debug("#%s mated #%s", id(self), id(other))
        return self, other

    def mutate(self, indpb):
        for gene in self.genes.values():
            gene.mutate(indpb)
        logger.debug("#%s mutated", id(self))
        return self,

    def similar(self, individual):
        logger.debug("Comparing RMSD between #%s and #%s",
                     id(self), id(individual))
        self.express()
        compound1 = next(gene for gene in self.genes.values()
                         if gene.__class__.__name__ == 'Molecule').compound
        atoms1 = sorted(compound1.mol.atoms, key=lambda x: x.serialNumber)
        coords1 = [a.coord() for a in atoms1]
        xf1 = compound1.mol.openState.xform
        self.unexpress()

        individual.express()
        compound2 = next(gene for gene in individual.genes.values()
                         if gene.__class__.__name__ == 'Molecule').compound
        atoms2 = sorted(compound2.mol.atoms, key=lambda x: x.serialNumber)
        coords2 = [a.coord() for a in atoms2]
        xf2 = compound2.mol.openState.xform
        individual.unexpress()

        sqdist = sum(xf1.apply(a).sqdistance(xf2.apply(b))
                     for a, b in zip(coords1, coords2))
        rmsd = math.sqrt(sqdist / ((len(coords1) + len(coords2)) / 2.0))
        logger.debug("RMSD: %f", rmsd)
        return rmsd < self.cfg.ga.similarity_rmsd

    def write(self, i):
        """
        # Maybe someday we can pickle it all :/
        filename = os.path.join(path, '{}_{}.pickle.gz'.format(name,i))
        with gzip.GzipFile(filename, 'wb') as f:
            cPickle.dump(self, f, 0)
        return filename
        """
        path = self.cfg.general.outputpath
        name = self.cfg.general.name
        COMPRESS = ZIP_DEFLATED if self.cfg.general.compress else ZIP_STORED
        self.express()
        zipfilename = os.path.join(path, '{}_{:03d}.zip'.format(name, i))
        with ZipFile(zipfilename, 'w', COMPRESS) as z:
            output = OrderedDict()
            for gene in self.genes.values():
                logger.debug("Writing %s to file", gene.name)
                filename = gene.write(path, "{}_{:03d}".format(name, i))
                if filename:
                    z.write(filename, os.path.basename(filename))
                    os.remove(filename)
                    output[gene.name] = os.path.basename(filename)
            try:
                output['score'] = list(self.fitness.values)
            except AttributeError:  # fitness not in individual :/
                raise
            z.writestr('{}_{:03d}.gaudi'.format(name, i),
                       yaml.dump(output, default_flow_style=False))
        self.unexpress()
        return zipfilename


class Fitness(deap.base.Fitness):

    """
    Augmented `Fitness` class to self-include `objectives` objects.

    It subclasses  DEAP's `Fitness` to include details of objectives being evaluated
    and a helper function to evaluate them all at once. Since Fitness it's an Attribute
    of every `individual`, it should result in a self-contained object.

    :parent:    A reference to the :class:`individual` that cointains the instance.
    :cache:     A reference to a mutable object that is maintained at `individual` level.
    :args:      Positional arguments that will be passed to `deap.base.Fitness.__init__`
    :kwargs:    Optional arguments that will be passed to `deap.base.Fitness.__init__`
    """

    objectives = OrderedDict()
    wvalues = ()

    def __init__(self, parent=None, cache=None, *args, **kwargs):
        self.parent = parent
        self.weights = self.parent.cfg.weights
        deap.base.Fitness.__init__(self, *args, **kwargs)
        self.cache = cache
        self.env = chimera.selection.ItemizedSelection()

        if not self.objectives:
            gaudi.plugin.load_plugins(self.parent.cfg.objectives,
                                      container=self.objectives,
                                      cache=self.cache, parent=self.parent,
                                      environment=self.env)

    def __deepcopy__(self, memo):
        copy_ = self.__class__(parent=self.parent)
        copy_.wvalues = self.wvalues + ()
        return copy_

    def evaluate(self):
        scores = []
        for name, obj in self.objectives.items():
            score = obj.evaluate()
            scores.append(score)
            logger.debug("%s fitness is %s", name, score)
        return scores
