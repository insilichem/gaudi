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

.. _developers:

================
Developers guide
================

Introduction to Genetic Algorithms
==================================

The GA in GaudiMM stands for Genetic Algorithm, a search heuristic inspired by natural selection that is used for optimization processes. 

Genetic Algorithms use a biologicist terminology. Each candidate solution to the problem is considered an **individual**, which is part of the so-called **population** (the set of all candidate solutions). 

The initial population is generated from scratch, almost always randomly. These individuals also comprise the first generation of the evolution process. As in nature, only the **fittest** survive. The survival process is simulated with an evaluation function, that tests them against the optimization variable(s). This is called **selection**. 

- How do we create an individual? With one or more **genes**, naturally. Genes describe how an individual should look like, of course!
- How do we evaluate that individual? With one or more **objectives**.

Also, as in nature, the fittest individuals are allowed to **mate** (exchange their defining values), and **mutate** (spontaneously modify their own defining values). This adds some more variability to the process, and that is key to survival.

After a number of generations repeating the process of selection, choosing the fittest over the weakest, we will obtain better and better solutions to our problem. Since it's all heuristics, it's up to us to stop at some point. We won't probably get *the* solution, but we can live with pretty good ones, right?

The One Max Problem
-------------------

Ok, that was a lot of biology and we are trying to code, I get it. Let's explain a classical GA example, the trivial **One Max Problem**, adapted from the `original deap documentation <https://deap.readthedocs.org/en/master/examples/ga_onemax.html>`_.

In this problem, we have a list of integers that can be either 0 or 1, and we want to obtain a list full of 1s. So, in this example, we have individuals defined by a single **gene** and evaluated with a single **objective**.

We build individuals with the gene **onemax**:

.. code-block:: python

    import random.randint
    def onemax(size=5):
        return [random.randint(0, 1) for i in range(size)]

.. code-block:: python
    
    adam = onemax(5) # returns [0, 0, 0, 1, 0]
    eve = onemax(5)  # returns [0, 1, 1, 0, 0]

The objective is also trivial. We have to maximize the sum of the numbers inside a given individual:

.. code-block:: python

    def evaluate(individual):
        return sum(individual)

So... which is one is fittest, ``adam`` or ``eve``? Obviously, ``eve``:

.. code-block:: python

    evaluate(adam)  # returns 1
    evaluate(eve)   # returns 2

Of course, an initial population is usually larger! At least, a hundred individuals. With such a trivial case, given a big enough population, we may obtain the solution in the first generation by pure change. However, we must not rely on the initial population as the only diversity source. 

Additional diversity is achieved with the mutation and mating operations, implemented as additional functions:

.. code-block:: python

    def mate(a, b):
        """ Let a and b mate, in hope of fitter children """
        i = crossover_point = random.random() * min(len(a), len(b))
        c, d = a[:], b[:]
        c[:i], d[i:] = d[:i], c[i:]
        return c, d

    def mutate(individual, probability):
        """ Spontaneous mutation at random places can result in a fitter individual """
        return [random.randint(0, 1) for i in individual if random.random() < probability]
        

Let's see how this is useful:

.. code-block:: python
    
    cain, abel = mate(adam, eve)
    # cain = [ 0, 1, 1, 1, 0 ]
    # abel = [ 0, 0, 0, 0, 0 ]
    evaluate(cain) # returns 3
    evaluate(abel) # returns 0


See? ``adam`` and ``eve`` gave birth to ``cain`` and ``abel``. ``cain`` had luck and inherited the good parts, while ``abel``... Well, he was not that lucky. In the next selection process, ``cain`` will be selected over ``abel``, and probably over its own father ``adam``. Now, the population (``cain`` and ``eve``) as a whole is fitter, with an average fitness of 2.5. That's higher than the average in the previous generation (1.5). Evolution!

Mutation works similarly:

.. code-block:: python

    enoch = mutate(cain)
    # enoch = [ 1, 1, 1, 1, 0]
    seth = mutate(eve)
    # seth = [ 0, 0, 1, 0, 0]

Take into account that mutations can be beneficial, like in the case of ``enoch``, but also detrimental, as in the case of ``seth``. Some of them will contribute to evolution, and some of them not. Lucky ones will be selected, the others, discarded.

By the way, deap already defines `some mutation and mating operators <https://deap.readthedocs.org/en/master/api/tools.html#operators>`_ for you that will work in most cases. So, hopefully, this part will be trivial.   

And that's it! Deap does the rest! So, to sum up, you only need to worry about:

- How to define your individuals.
- How to evaluate them.
- How to implement mutation and mating (normally, with deap built-in operators).

If you want to know more about Deap and Genetic Algorithms, go check their `documentation <https://deap.readthedocs.org/en/master/index.html>`_. It's great!

Our implementation
==================

GaudiMM is built as an extensible and highly modular Python platform. Although the main focus is Chemistry and molecular design, you can use your own genes and objectives. You can think of GaudiMM as a new API for `deap <https://github.com/deap/deap>`_ that provides an object-oriented interface to easily create new individuals and objectives.

In ``deap`` an individual can be any Python object (check their `overview <https://deap.readthedocs.org/en/master/overview.html>`_ and `GA examples <https://deap.readthedocs.org/en/master/examples/ga_onemax.html>`_), which is a very versatile approach, but it tends to be very limited when your individual gets complex. For example, if an individual needs to be defined by several genes with different mutation strategies.

In GaudiMM, each **individual** is a :class:`gaudi.base.Individual`, which is a very (bio)fancy name for a list of ``genes``. To create a ``gene``, you just subclass :class:`gaudi.genes.GeneProvider` and define the needed methods: ``express``, ``unexpress``, ``mutate``, and ``mate``. The :class:`gaudi.base.Individual` class then provides some wrapper methods that call the respective counterparts in each ``gene``.

To evaluate the fitness of an individual, you must first define the set of evaluation functions. Each function is called ``objective``, and you keep them inside a :class:`gaudi.base.Environment`.

To create a new ``objective``, you have to subclass :class:`gaudi.objectives.ObjectiveProvider`, which provides a very simple interface: ``evaluate``. Define your function there, and that's it!

.. todo::

    * Tutorial: How to create your own gene
    * Tutorial: How to create your own objective

