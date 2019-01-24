%PYTHON% -m pip install prody==1.8.2
if errorlevel 1 exit 1
%PYTHON% -m pip install . --no-deps -v
