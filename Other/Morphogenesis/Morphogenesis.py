#!/usr/bin/env python
##############################################################################
# EVOLIFE  www.dessalles.fr/Evolife                    Jean-Louis Dessalles  #
#            Telecom ParisTech  2014                       www.dessalles.fr  #
##############################################################################

##############################################################################
# Basic implementation of Turing's reaction-diffusion model                  #
##############################################################################


""" Basic implementation of Turing's reaction-diffusion model              
"""

import sys
sys.path.append('..')
sys.path.append('../../..')

import Evolife.Ecology.Observer				as EO
import Evolife.Scenarii.Parameters 			as EPar
import Evolife.QtGraphics.Evolife_Window 	as EW
import Landscapes

from Evolife.Tools.Tools import boost
print boost()   # A technical trick that sometimes provides impressive speeding up

	
import random
import math


# global functions
Quantify = lambda x, Step: int(x / Step) * Step if Step else x

def Limitate(x, Min=0, Max=1, Step=0): 
	if Step == 0:	Step = Gbl.P('Precision')
	return min(max(Quantify(x, Step),Min), Max)

def BlueShade(C, Min=0, Max=1):	
	C1 = Limitate(C, Min=Min, Max=Max)
	# Blue colours between 34 and 45
	# return 45 - int((10 * C1/Max))		# 10 blue shades
	if int(10*C1):	return 45 - int((10 * C1/Max))		# 10 blue shades
	# if int(3*C1):	return 45 - int((3 * C1))		# 3 blue shades
	else:	return 2	# white colour == invisible

def RedShade(C, Min=0, Max=1):	
	# Red colours between 22 and 33
	C1 = Limitate(C, Min=Min, Max=Max)
	if int(8*C1):	return 33 - int(8 * C1/Max)	# 8 red shades
	else:	return 2	# white colour == invisible



class Scenario(EPar.Parameters):
	def __init__(self):
		# Parameter values
		EPar.Parameters.__init__(self, CfgFile='_Params.evo')	# loads parameters from configuration file

	def P(self, ParamName):	return self.Parameter(ParamName)
		
class Morphogenesis_Observer(EO.Observer):
	""" Stores global variables for observation
	"""
	def __init__(self, Scenario):
		EO.Observer.__init__(self, Scenario) # stores global information
		self.CurrentChanges = []	# stores temporary changes

	def record(self, Info):
		# stores current changes
		# Info is a couple (InfoName, Position) and Position == (x,y)
		self.CurrentChanges.append(Info)

	def get_data(self, Slot):
		" this is called when display is required "
		if Slot == 'Positions':
			CC = self.CurrentChanges
			self.CurrentChanges = []
			return tuple(CC)
		else:	return EO.Observer.get_data(self, Slot)
			
			
		
class LandCell(Landscapes.LandCell):
	" Defines what is stored at a given location "

	def __init__(self):
		# Cell content is defined as a couple  (ConcentrationA, ConcentrationB)
		Landscapes.LandCell.__init__(self, (0, 0), VoidCell=(0, 0))
		self.Colour = 2	# white == invisible
		self.OldColour = None

	def activated(self, Future=False):
		" tells whether a cell is active "
		return self.Content(Future=Future)[1] > 0
	
	def Update(self):
		self.Present = self.Future	# erases history
		return not self.activated()
		
	def getColour(self, Max=1):
		" return color corresponding to content, but only if it has changed "
		self.Colour = BlueShade(self.Content()[1], Max=Max)
		Col = None
		if self.Colour != self.OldColour:	Col = self.Colour
		self.OldColour = self.Colour
		return Col
	
	#------------------------#
	# Reaction               #
	#------------------------#
	# Gray-Scott model       #
	#------------------------#
	def Reaction(self, F, k):
		"	Reaction between local products A and B "
		(Ca0, Cb0) = self.Content()
		deltaB = Ca0 * Cb0 * Cb0
		deltaA = F * (1 - Ca0) - deltaB
		deltaB -= (F + k) * Cb0
		return ((Limitate(Ca0 + deltaA, Max=Gbl.P('MaxA'))), Limitate(Cb0 + deltaB, Max=Gbl.P('MaxB')))
		
	
class Landscape(Landscapes.Landscape):
	" Defines a 2D square grid "

	def __init__(self, Size, DiffusionCoefficients, NeighbourhoodRadius):
		Landscapes.Landscape.__init__(self, Width=Size, CellType=LandCell)	# calls local LandCell definition
		# Computing actual diffusion coefficients
		self.NeighbourhoodRadius = NeighbourhoodRadius
		self.DC = DiffusionCoefficients

	def Seed(self, Center, Value, Radius=5):
		" Creation of a blob "
		for Pos1 in self.neighbours(Center, Radius):	# immediate neighbours
			self.Modify(Pos1, Value, Future=False)
			
	#------------------------#
	# Diffusion              #
	#------------------------#
	def activate(self, Pos0):
		" Cell located at position 'Pos0' produces its effect on neighbouring cells "
		# This function is called by 'activate'
		(Ca0, Cb0) = self.Cell(Pos0).Content()	# actual concentration values
		Neighbours = self.neighbours(Pos0, self.NeighbourhoodRadius)
		##### First method: each cell pours part of its content into its neighbours
		(depletionA, depletionB) = (0,0)
		for Pos1 in Neighbours:	# immediate neighbours
			(Ca1, Cb1) = self.Content(Pos1, Future=True)	# future concentration values
			deltaA = min(self.DC[0] * Ca0/len(Neighbours), 1 - Ca1)	# saturation
			deltaB = min(self.DC[1] * Cb0/len(Neighbours), 1 - Cb1)	# saturation
			self.Modify(Pos1, (	Limitate(Ca1 + deltaA, Max=Gbl.P('MaxA')),
								Limitate(Cb1 + deltaB, Max=Gbl.P('MaxB'))),
								Future=True)	# add new activation 
			depletionA += deltaA
			depletionB += deltaB
		self.Modify(Pos0, (	Limitate(Ca0 - depletionA, Max=Gbl.P('MaxA')), 
							Limitate(Cb0 - depletionB, Max=Gbl.P('MaxB'))),
							Future=True)
		##### Second method: averaging local gradients - much faster, but probably inaccurate
		# # (Ca0f, Cb0f) = self.Cell(Pos0).Content(Future=True)	# future concentration values
		# # (AvgA, AvgB) = (0, 0)
		# # for Pos1 in Neighbours:	# neighbours
			# # (Ca1, Cb1) = self.Content(Pos1)	# old concentration values
			# # AvgA += Ca0 - Ca1
			# # AvgB += Cb0 - Cb1
		# # AvgA = AvgA * self.DC[0] / self.NeighbourhoodSize
		# # AvgB = AvgB * self.DC[1] / self.NeighbourhoodSize
		# # self.Modify(Pos0, (Limitate(Ca0f - AvgA, Max=Gbl.P('MaxA')), Limitate(Cb0f - AvgB, Max=Gbl.P('MaxB'))),
					# # Future=True)

def One_Step():
	""" This function is repeatedly called by the simulation thread.
		One agent is randomly chosen and decides what it does
	"""

	Observer.season()	# increments StepId
	dotSize = Gbl.P('DotSize')
	if Observer.Visible():	
		# for (Position, Cell) in Land.travel():
		for Position in Land.ActiveCells.list():
			Cell = Land.Cell(Position)
			# displaying concentrations
			Colour = Cell.getColour(Max=Gbl.P('MaxB'))
			if Colour is not None:
				Observer.record(('C%d_%d' % Position, Position + (Colour, dotSize))) 

	#############
	# Diffusion #
	#############
	Land.activation()	# diffusion
	Land.update()	# Let's update concentrations
	
	#############
	# Reaction  #
	#############
	for (Position, Cell) in Land.travel():
		Land.Modify(Position, Cell.Reaction(Gbl.P('F'), Gbl.P('k')), Future=False)

	if len(Land.ActiveCells) == 0: return False

		  
	return True
			


if __name__ == "__main__":
	print __doc__

	
	#############################
	# Global objects			#
	#############################
	Gbl = Scenario()
	Observer = Morphogenesis_Observer(Gbl)	  # Observer contains statistics
	Land = Landscape(Gbl.P('LandSize'), (Gbl.P('Da'), Gbl.P('Db')), Gbl.P('NeighbourhoodRadius'))	  # 2D square grid
	Land.Seed((Gbl.P('LandSize')//3, Gbl.P('LandSize')//3), (0.4, 0.9), Radius=4)
	Land.Seed((2*Gbl.P('LandSize')//3, 2*Gbl.P('LandSize')//3), (0.4, 0.9), Radius=4)
	# print len(Land.ActiveCells)
	
	# Land.setAdmissible(range(101))	# max concentration
	# Observer.recordInfo('BackGround', 'white')
	Observer.recordInfo('FieldWallpaper', 'white')
	
	EW.Start(One_Step, Observer, Capabilities='RP')

	print "Bye......."
	
__author__ = 'Dessalles'
