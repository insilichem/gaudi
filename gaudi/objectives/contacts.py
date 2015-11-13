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
This objective provides a wrapper around Chimera's
`DetectClash` that detects clashes and contacts.

Clashes are understood as steric conflicts that increases the energy
of the system. They are evaluated as the sum of volumetric overlapping
of the Van der Waals' spheres of the implied atoms.

Contacts are considered as stabilizing, and they are evaluated with a
Lennard-Jones 12-6 like function.

"""

# Python
import logging
# Chimera
import chimera
import DetectClash
import ChemGroup as cg
# GAUDI
from gaudi.objectives import ObjectiveProvider

logger = logging.getLogger(__name__)


def enable(**kwargs):
    return Contacts(**kwargs)


class Contacts(ObjectiveProvider):

    """
    Contacts class

    Parameters
    ----------
    probe : str
        Name of molecule gene that is object of contacts analysis
    radius : float
        Maximum distance from any point of probe that is searched
        for possible interactions
    which : {'hydrophobic', 'clashes'}
        Type of interactions to measure
    threshold : float, optional
        Maximum overlap of van-der-Waals spheres that is considered as
        a contact (attractive). If the overlap is greater, it's 
        considered a clash (repulsive)
    threshold_h : float, optional
        If the involved atoms can form a H-bond, use this threshold instead.
    treshold_c : float, optional
        If the involved atoms can form a hydrophobic patch, use this threshold 
        instead. (sure??)
    cutoff : float, optional
        If the overlap volume is greater than this, a penalty is applied. 
        Useful to filter bad solutions.
    """

    def __init__(self, probe=None, radius=5.0, which='hydrophobic',
                 threshold=0.6, threshold_h=0.2, threshold_c=0.6, cutoff=100.0,
                 *args, **kwargs):
        ObjectiveProvider.__init__(self, **kwargs)
        self.which = which
        self.radius = radius
        self.threshold = threshold
        self.threshold_h = threshold_h
        self.threshold_c = threshold_c
        self.cutoff = cutoff
        self._probe = probe

    def molecules(self, ind):
        return tuple(m.compound.mol for m in ind.genes.values()
                     if m.__class__.__name__ == "Molecule")

    def probe(self, ind):
        return ind.genes[self._probe].compound.mol

    def evaluate(self, ind):
        """
        Get interacting pairs of atoms with DetectClash (Chimera's) and submit them
        to our parsers.
        """
        test_atoms = self._surrounding_atoms(ind)
        clashes = DetectClash.detectClash(test_atoms, test=test_atoms, intraRes=True,
                                          interSubmodel=True, clashThreshold=self.threshold,
                                          hbondAllowance=self.threshold_h, assumedMaxVdw=2.1,
                                          bondSeparation=4)
        try:
            positive, negative = self._parse_clashes_c(clashes, ind)
        except StopIteration:
            positive, negative = 0, 0

        if self.which == 'clashes':
            clashscore = sum(abs(a[3]) for a in negative) / 2
            if clashscore > self.cutoff:
                clashscore = -1000 * self.weight
            return clashscore
        elif self.which == 'hydrophobic':
            return sum(1 - a[3] for a in positive) / 2

    ###
    def _parse_clashes(self, clashes, ind):
        ALIPH = ['C3', [[cg.C, [cg.C, [cg.C, [cg.R, cg.R, cg.R, cg.R]],
                                cg.R, cg.R]], cg.R, cg.R, cg.R]], [1, 1, 1, 1, 1, 0, 0]
        AROMATIC = set(
            a for g in cg.findGroup("aromatic ring", self.molecules(ind)) for a in g)
        ALIPHATIC = set(a for g in cg.findGroup(ALIPH, self.molecules(ind))
                        for a in g if a not in AROMATIC)

        positive, negative = [], []
        for a1, c in clashes.items():
            for a2, dist in c.items():
                if dist <= self.threshold and a1.residue != a2.residue:
                    if a1 in AROMATIC and a2 in AROMATIC:
                        positive.append(
                            [a1, a2, dist, self._lennard_jones(a1, a2)])
                    elif a1 in AROMATIC and a2 in ALIPHATIC:
                        positive.append(
                            [a1, a2, dist, self._lennard_jones(a1, a2)])
                    elif a1 in ALIPHATIC and a2 in AROMATIC:
                        positive.append(
                            [a1, a2, dist, self._lennard_jones(a1, a2)])
                elif dist > self.threshold:
                    negative.append(
                        [a1, a2, dist, self._vdw_vol_overlap(a1, a2)])

        return positive, negative

    def _parse_clashes_c(self, clashes, ind):
        """
        Interpret contacts provided by DetectClash.

        Parameters
        ----------
        clashes : dict of dict
            Output of DetectClash. It's a dict of atoms, whose values are dicts.
            These subdictionaries contain all the contacting atoms as keys, and
            the respective distance as values.
        ind : Individual

        Returns
        -------
        positive : list of list
            Each sublist depict an interaction, with four items: the two involved
            atoms, their distance, and their Lennard-Jones score.
        negative : list of list
            Each sublist depict an interaction, with four items: the two involved
            atoms, their distance, and their volumetric overlap.

        .. note ::
            First, collect atoms that can be involved in hydrophobic interactions.
            Namely, C and S.

            Then, iterate the contacting atoms, getting the distances. For each
            interaction, analyze the distance and, based on the threshold, determine
            if it's attractive or repulsive.

            Attractive interactions are weighted with a Lennard-Jones like function
            (``_lennard_jones``), while repulsive attractions are measured with
            the volumetric overlap of the involved atoms' Van der Waals spheres.

        """
        vdwatoms = set(
            a for m in self.molecules(ind) for a in m.atoms if a.element.name in ('C', 'S'))

        positive, negative = [], []
        for a1, c in clashes.items():
            for a2, dist in c.items():
                if dist <= self.threshold and a1.molecule != a2.molecule:
                    if a1 in vdwatoms and a2 in vdwatoms:
                        positive.append(
                            [a1, a2, dist, self._lennard_jones(a1, a2)])
                elif dist > self.threshold:
                    negative.append(
                        [a1, a2, dist, self._vdw_vol_overlap(a1, a2)])

        return positive, negative

    def _surrounding_atoms(self, ind):
        """
        Get atoms in the search zone, based on the molecule and the radius
        """
        self.zone.clear()
        self.zone.add(self.probe(ind).atoms)
        self.zone.merge(chimera.selection.REPLACE,
                        chimera.specifier.zone(
                            self.zone, 'atom', None, self.radius, self.molecules(
                                ind)
                        )
                        )
        return self.zone.atoms()

    @staticmethod
    def _lennard_jones(a1, a2):
        """
        VERY rough aproximation of a Lennard-Jones score (12-6).

        Parameters
        ----------
        a1, a2 : chimera.Atom

        .. todo ::

            The distance is usually computed in the clash parsers, so get it
            instead of computing it again. At least, as an optional kw.
        """
        dist = a1.xformCoord().distance(a2.xformCoord())
        zero = 0.98 * (a1.radius + a2.radius)
        x = zero / dist
        return (x ** 12 - 2 * x ** 6)

    @staticmethod
    def _vdw_vol_overlap(a1, a2):
        """
        Volumetric overlap of Van der Waals spheres of atoms.

        Parameters
        ----------
        a1, a2 : chimera.Atom

        .. note ::
            Adapted from Eran Eyal, Comput Chem 25: 712-724, 2004
        """
        PI = 3.14159265359
        d = a1.xformCoord().distance(a2.xformCoord())
        if not d:
            return 1000
        h_a, h_b = 0, 0
        if d and d < a1.radius + a2.radius:
            h_a = (a2.radius ** 2 - (d - a1.radius) ** 2) / (2 * d)
            h_b = (a1.radius ** 2 - (d - a2.radius) ** 2) / (2 * d)
        v = (PI / 3) * (h_a ** 2) * (3 * a1.radius - h_a) + \
            (PI / 3) * (h_b ** 2) * (3 * a2.radius - h_b)
        return v
