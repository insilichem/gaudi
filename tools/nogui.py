#!chimera

## VERSION 2
# Bat and Sh files are not needed anymore
# Everything is handled through Python

#import chimera
from chimera import runCommand as rc
from chimera import UserError, pathFinder, dialogs, replyobj
import os
from sys import argv
from time import strftime, sleep
import subprocess, shlex

# Run NO GUI Chimera
if len(argv) == 1:
	raise UserError("ERROR! Please provide the script and arguments")

timestamp = strftime("%Y%m%d%H%M%S")

wd = os.path.dirname(os.path.realpath(argv[0])).replace('\\', '/') + '/'
chimera_path = pathFinder().dataRoot.replace('share','bin')
exe = '.exe' if os.name == 'nt' else ''
session = wd + "tmpses/session" + timestamp
script = ' '.join(argv[1:]) # 1st argument is this script path

# Close all dialogs before saving session
# for d in dialogs.activeDialogs():
# 	d.Close()

# Save current session with tmp file
if not os.path.isdir(wd + "tmpses/"):
	os.makedirs(wd + "tmpses/")
rc("save {0}".format(session))

command = '{0}/chimera{1} --debug --nogui --script "{2}" "{3}.py"'.format(chimera_path,exe,script,session)
print "Run this command in terminal if you are having problems with your script:\n" + command

run = subprocess.Popen(shlex.split(command), stdout=subprocess.PIPE)
replyobj.status("Script '" + script + "' now running.", blankAfter=1)
#run.wait() #Uncomment if progress is not in use
# Progress reporting
for line in iter(run.stdout.readline, ""):
	if 'Progress' in line:
		replyobj.status("Script: " + str(line))
replyobj.status("Script '" + script + "' was executed.", blankAfter=5)

