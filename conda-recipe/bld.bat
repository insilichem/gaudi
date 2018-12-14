%PYTHON% -m pip install pychimera prody==1.8.2 -v
if errorlevel 1 exit 1
%PYTHON% -m pip install . --no-deps -v
