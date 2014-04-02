@ECHO off
ECHO This will run a script in Chimera in debug mode (no GUI)
SET chimera=C:\Program Files\Chimera 2014-01-10

ECHO Running script and args: %1
ECHO in session %2

cd "%chimera%"
"%chimera%\bin\chimera.exe" --debug --nogui --script %1  %2