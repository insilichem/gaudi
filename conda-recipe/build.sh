#!/bin/bash
$PYTHON -m pip install pychimera munch voluptuous click boltons deap pyyaml prody==1.8.2 -vv
$PYTHON -m pip install . --no-deps -vv