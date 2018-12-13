%PYTHON% -m pip install pychimera munch voluptuous click boltons deap pyyaml prody==1.8.2 -vv
if errorlevel 1 exit 1
%PYTHON% -m pip install . --no-deps -vv
