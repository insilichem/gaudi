Developers guide
================

GAUDI is built as an extensible and highly modular Python platform. Although the main focus is Chemistry and molecular design, you can use your own genes and objectives. You can think of GAUDI as a new API for `deap <https://github.com/deap/deap>`_ that provides an object-oriented interface to easily create new individuals and objectives.


In ``deap`` an individual can be any Python object, which is a very versatile approach, but it tends to be very limited when your individual gets complex: several different genes with different requirements. What if I want to use a mutation operator on some genes of the individual, and a different mutation system on another gene?

In GAUDI, each **individual** is a :class:`gaudi.base.Individual`, which is a very (bio)fancy name for a set of ``genes``. To create a ``gene``, you just subclass :class:`gaudi.genes.GeneProvider` and define the needed methods: express and unexpress, mutate, and mate. The :class:`gaudi.base.Individual` class then provides some wrapper methods that call the respective counterparts in each ``gene``.

To evaluate the fitness of an individual, you must first define the set of evaluation functions. Each function is called ``objective``, and you keep them inside a :class:`gaudi.base.Environment`.

To create a new ``objective``, you have to subclass :class:`gaudi.objectives.ObjectiveProvider`, which provides a very simple interface: ``evaluate``. Define your function there, and that's it!

.. todo::

    * Tutorial: How to create your own gene
    * Tutorial: How to create your own objective
