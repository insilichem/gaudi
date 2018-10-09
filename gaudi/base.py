#!/usr/bin/env python
# -*- coding: utf-8 -*-

##############
# GaudiMM: Genetic Algorithms with Unrestricted
# Descriptors for Intuitive Molecular Modeling
#
# https://github.com/insilichem/gaudi
#
# Copyright 2017 Jaime Rodriguez-Guerra, Jean-Didier Marechal
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
##############

"""
Contains the core classes we use to build individuals
(potential solutions of the optimization process).

"""

# Python
from collections import OrderedDict
from contextlib import contextmanager
from copy import deepcopy
from zipfile import ZipFile, ZIP_DEFLATED, ZIP_STORED
import logging
import os
import pprint
import sys
# Chimera
import chimera
# External dependencies
import deap.base
import yaml
# GAUDI
import gaudi.plugin
import gaudi.similarity
from gaudi.exceptions import MoleculesNotFound

pp = pprint.PrettyPrinter(4)
logger = logging.getLogger(__name__)


class BaseIndividual(object):

    """
    Base class for `individual` objects that are evaluated by DEAP.

    Each individual is a potential solution. It contains all that is needed
    for an evaluation. With multiprocessing in mind, individuals should be
    self-contained so it can be passed between threads.

    The defined methods are only wrapper calls to the respective methods of each
    gene.

    Parameters
    ----------
    cfg : gaudi.parse.Settings
        The full parsed object from the configuration YAML file.
    cache : dict or dict-like
        A mutable object that can be used to store values across instances.
    dummy : bool
        If True, create an uninitialized Individual, only containing the
        cfg attribute. If false, call `__ready__` and complete initialization.

    Attributes
    ----------
    __CACHE : dict
        Class attribute that caches gene data across instances
    __CACHE_OBJ : dict
        Class attribute that caches objectives data across instances


    .. todo::

        :meth:`write` should use Pickle and just save the whole object, but
        Chimera's inmutable objects (Atoms, Residues, etc) get in the way. A
        workaround may be found if we take a look a the session saving code.

    """

    _CACHE = {}
    _CACHE_OBJ = {}

    def __init__(self, cfg=None, cache=None, **kwargs):
        logger.debug("Creating new individual with id %s", id(self))
        self.cfg = cfg
        self.genes = OrderedDict()
        self.fitness = None
        self.expressed = False
        self.__ready__()
        self.__expression_hooks__()

    def __ready__(self):
        """
        A *post-init* method used to avoid initialization problems with
        `__deepcopy__`. It's just the second part of a three-stage
        `__init__`.
        """
        if self.cfg is not None:
            gaudi.plugin.load_plugins(self.cfg.genes, container=self.genes,
                                      parent=self,
                                      cx_eta=self.cfg.ga.cx_eta,
                                      mut_eta=self.cfg.ga.mut_eta,
                                      mut_indpb=self.cfg.ga.mut_indpb)

            self.fitness = Fitness(self.cfg.weights)
            mod, fn = self.cfg.similarity.module.rsplit('.', 1)
            self._similarity = getattr(sys.modules[mod], fn)

    def __expression_hooks__(self):
        """
        Third part of initialization. This method is run after expression
        of almost-built individuals to retrieve information that would
        only be available at expression time.
        """
        with expressed(self):
            for name, gene in self.genes.items():
                gene.__expression_hooks__()

    def __deepcopy__(self, memo):
        new = self.__class__(cfg=self.cfg)
        new.genes = deepcopy(self.genes, memo)
        new.fitness = deepcopy(self.fitness, memo)
        new._similarity = self._similarity
        new.expressed = self.expressed
        return new

    def evaluate(self, environment):
        """
        Express individual, evaluate it and unexpress it.

        Parameters
        ----------
        environment : Environment
            Objectives that will evaluate the individual
        """
        logger.debug("Evaluating individual #%s", id(self))
        self.express()
        score = environment.evaluate(self)
        self.unexpress()
        return score

    def pre_express(self):
        pass

    def express(self):
        """
        Express genes in this environment. Very much like 'compiling' the
        individual to a chimera.Molecule.
        """
        self.pre_express()
        for name, gene in self.genes.items():
            logger.debug("Expressing gene %s with allele\n%s",
                         name, pp.pformat(gene.allele))
            gene.express()
        self.post_express()
        self.expressed = True

    def post_express(self):
        pass

    def pre_unexpress(self):
        pass

    def unexpress(self):
        """
        Undo .express()
        """
        self.pre_unexpress()
        for gene in reversed(self.genes.values()):
            logger.debug("Reverting expression of gene %s", gene.name)
            gene.unexpress()
        self.post_unexpress()
        self.expressed = False

    def post_unexpress(self):
        pass

    def mate(self, other):
        """
        Recombine genes of `self` with `other`. It simply calls `mate` on
        each gene instance

        Parameters
        ----------
        other : Individual
            Another individual to mate with.
        """
        for this_gene, other_gene in zip(self.genes.values(), other.genes.values()):
            this_gene.mate(other_gene)
        logger.debug("#%s mated #%s", id(self), id(other))
        return self, other

    def mutate(self, indpb):
        """
        Trigger a round of possible mutations across all genes

        Parameters
        ----------
        indpb : float
            Probability of suffering a mutation
        """
        for gene in self.genes.values():
            gene.mutate(indpb)
        logger.debug("#%s mutated", id(self))
        return self,

    def similar(self, other):
        """
        Compare `self` and `other` with a similarity function.

        Returns
        -------
        bool
        """
        return self._similarity(self, other,
                                *self.cfg.similarity.args,
                                **self.cfg.similarity.kwargs)

    def write(self, i, path=None):
        """
        Export the individual to a mol2 file

        Parameters
        ----------
        i : int
            Individual identificator in current generation or hall of fame

        .. note ::

            Maybe someday we can pickle it all :/
            >>> filename = os.path.join(path, '{}_{}.pickle.gz'.format(name,i))
            >>> with gzip.GzipFile(filename, 'wb') as f:
            >>>     cPickle.dump(self, f, 0)
            >>> return filename
        """
        path = path if path is not None else self.cfg.output.path
        name = self.cfg.output.name
        COMPRESS = ZIP_DEFLATED if self.cfg.output.compress else ZIP_STORED
        self.express()
        zipfilename = os.path.join(path, '{}_{:03d}.zip'.format(name, i))
        with ZipFile(zipfilename, 'w', COMPRESS) as z:
            output = OrderedDict()
            for name, gene in self.genes.items():
                logger.debug("Writing %s to file", name)
                filename = gene.write(path, "{}_{:03d}".format(name, i))
                if filename:
                    z.write(filename, os.path.basename(filename))
                    os.remove(filename)
                    output[name] = os.path.basename(filename)
            output['score'] = list(self.fitness.values)
            z.writestr('{}_{:03d}.gaudi'.format(name, i),
                       yaml.dump(output, default_flow_style=False))
        self.unexpress()
        return zipfilename

    def clear_cache(self):
        self._CACHE.clear()
        self._CACHE_OBJ.clear()
        for gene in self.genes.values():
            gene.clear_cache()


class MolecularIndividual(BaseIndividual):


    def __init__(self, *args, **kwargs):
        self._molecules = OrderedDict()
        BaseIndividual.__init__(self, *args, **kwargs)

    def __ready__(self):
        BaseIndividual.__ready__(self)
        for name, gene in self.genes.items():
            gene.__ready__()
            if gene.__class__.__name__ == 'Molecule':
                self._molecules[name] = gene

    def __deepcopy__(self, memo):
        new = self.__class__(cfg=self.cfg)
        new.genes = deepcopy(self.genes, memo)
        new.fitness = deepcopy(self.fitness, memo)
        new._similarity = self._similarity
        new.expressed = self.expressed
        for name, gene in new.genes.items():
            gene.parent = new
            if gene.__class__.__name__ == 'Molecule':
                new._molecules[name] = gene
        return new

    def post_express(self):
        for m in self._molecules.values():
            m._expressed_coordinates = m.xyz()

    def find_molecule(self, name):
        try:
            return self._molecules[name]
        except KeyError:
            raise MoleculesNotFound('Molecule {} not found. '
                                    'Try with {}'.format(name, sorted(self._molecules.keys())))
    def xyz(self, gene=None):
        try:
            return self._molecules[gene]._expressed_coordinates
        except KeyError:
            return {name: mol._expressed_coordinates
                    for name, mol in self._molecules.items()}


@contextmanager
def expressed(*individuals):
    try:
        for individual in individuals:
            individual.express()
        yield individuals
    finally:
        for individual in individuals:
            individual.unexpress()


class Environment(object):

    """
    Objective container and helper to evaluate an individual. It must be
    instantiated with a gaudi.parse.Settings object.

    Parameters
    ----------
    cfg : gaudi.parse.Settings
        The parsed configuration YAML file that contains objectives information
    """

    def __init__(self, cfg=None, *args, **kwargs):
        self.zone = chimera.selection.ItemizedSelection()
        self.objectives = OrderedDict()
        self.cfg = cfg
        if self.cfg is not None:
            self.weights = self.cfg.weights
            gaudi.plugin.load_plugins(self.cfg.objectives,
                                      container=self.objectives,
                                      zone=self.zone, environment=self,
                                      precision=self.cfg.output.precision)

    def evaluate(self, individual):
        """
        individual : Individual
        """
        scores = []
        with expressed(individual):
            for name, obj in self.objectives.items():
                score = obj.evaluate(individual)
                score = round(score, obj.precision)
                scores.append(score)
                logger.debug("%s fitness is %s", name, score)
        return scores

    def clear_cache(self):
        for obj in self.objectives.values():
            obj.clear_cache()


class Fitness(deap.base.Fitness):

    wvalues = ()

    def __init__(self, weights):
        self.weights = weights
        deap.base.Fitness.__init__(self)

    def __deepcopy__(self, memo):
        new = self.__class__(self.weights)
        new.wvalues = self.wvalues + ()
        return new
