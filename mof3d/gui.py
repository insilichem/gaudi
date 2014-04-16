#!/usr/bin/python

# MOF3D
# Multi-Objective Force-Field-Free Docking
# A docking module for UCSF Chimera
# Jaime RGP <https://bitbucket.org/jrgp> @ UAB, 2014

import chimera
from chimera.baseDialog import ModelessDialog
import Tkinter, Pmw

class Mof3dDialog(ModelessDialog):
	name = "mof3d"
	title = "MOF3D / Multi-Objective Force-field free Docking"
	provideStatus = True
	statusPosition = "above"
	help = ""

	def fillInUI(self, parent):
		pass

	def OK(self, *args):
		ModelessDialog.OK(self,*args)
	def Apply(self, *args):
		pass
	def Close(self):
		ModelessDialog.Close(self)
