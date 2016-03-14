#!/usr/bin/env python
##############################################################################
# EVOLIFE  www.dessalles.fr/Evolife                    Jean-Louis Dessalles  #
#            Telecom ParisTech  2014                       www.dessalles.fr  #
##############################################################################

##############################################################################
# Example to show how to use Evolife's graphic system                        #
##############################################################################

""" This example shows how to use Evolife's graphic system (based on PyQT)
	to run simulations that display images,dots and lines, and curves.
"""

#	-----------------------------------------------------------------------------------
#	Evolife provides a window system that you can use for you own simulation.
#
#	* Simulation:
#	-	Re-implement the function 'One_Run' in the class 'Population'.
#		This function is repeatedly called by Evolife. It should perform the simulation.
#	-	This function typically executes individual behaviour: in this example, 
#		function 'move'. You may just re-write the latter.
#
#	* Display
#	Evolife gets instructions for display through the class 'Observer'.
#	Observer should return appropriate data to Evolife's requests, as indicated.
#	This means that the simulation must keep Observer informed of relevant changes.
#
#	* Starting Evolife
#		Evolife can be started with various capabilities for display 
#		(curves, dots and lines, genomes, links...) as indicated.
#	-----------------------------------------------------------------------------------

from time import sleep
from random import randint

import sys
sys.path.append('../..')

import Evolife.QtGraphics.Evolife_Window as EW
import Evolife.Ecology.Observer as EO



# Global constants
class Global:
	def __init__(self):
		self.NbAgents = 100
		self.TimeLimit = 1000

Gbl = Global()


class Observer(EO.Observer):
	""" Stores all values that should be displayed
		May also store general information
	"""

	def __init__(self):
		EO.Observer.__init__(self)
		self.TimeLimit = Gbl.TimeLimit
		self.Positions = dict() # position of agents, for display
		
	def get_data(self, Slot):

		# This function is called each time the window wants to update display
		
		#--------------------------------------------------------------------------#
		# Moving agents or drawing lines	(option F or R when starting Evolife)  #
		#--------------------------------------------------------------------------#
		if Slot == 'Positions':	
			# Should return a tuple or list of coordinates
			# Two formats are recognized:
			#	- ((Agent1Id, Coordinates1), (Agent2Id, Coordinates2), ...)
			#	- (Coordinates1, Coordinates2, ...)
			# The first format allows to move and erase agents
			# The second format is merely used to draw blobs and lines
			# Coordinates have the following form:
			#	(x, y, colour, size, ToX, ToY, segmentColour, thickness)
			#	(shorter tuples are automatically continued with default values)
			# The effect is that a blob of size 'size' is drawn at location (x,y) (your coordinates, not pixels)
			# and a segment starting from that blob is drawn to (ToX, ToY) (if these values are given)
			# If you change the coordinates of an agent in the next call, it will be moved.
			#
			# There are two modes, 'Field' and 'Region' (see below option F and R at Evolife's start)
			# 	- In the 'Field' mode, all agents should be given positions at each call.
			# 	  Missing agents are destroyed from display.
			#	- In the 'Region' mode, you may indicates positions only for relevant agents
			#	  To destroy an agent from display, give a negative value to its colour.
		
			if self.StepId == 0:
				# drawing a green square
				return ((5, 5, 'green', 0, 95, 5, 'green', 3),
						(95, 5, 'green', 8, 95, 95, 'green', 3),
						(95, 95, 'green', 0, 5, 95, 'green', 3),
						(5, 95, 'green', 8, 5, 5, 'green', 3))
			return self.Positions.items()
		#--------------------------------------------------------------------------#
		# Displaying genomes	(option G when starting Evolife)                   #
		#--------------------------------------------------------------------------#
		elif Slot == 'DNA':
			# Should return a list (or tuple) of genomes
			# Each genome is a tuple like (0,1,1,0,1,...)
			# All genomes should have same length
			return None
			
		#--------------------------------------------------------------------------#
		# Displaying social links	(option N (== network) when starting Evolife)  #
		#--------------------------------------------------------------------------#
		elif Slot == 'Network':
			# Should return a list (or tuple) of friends
			# 	((Agent1Id, Friends1), (Agent2Id, Friends2), ...)
			# Agents' Ids should be consistent with 'Positions'
			# Friends = list of agents Ids to which the agent is connected
			# Currently, only links to best friends are displayed
			return None
			
		else:   return EO.Generic_Observer.get_data(self, Slot)	# basic behaviour; 

	def get_info(self, Slot):
		#--------------------------------------------------------------------------#
		# Trajectories (option T when starting Evolife)                            #
		#--------------------------------------------------------------------------#
		if Slot == 'Trajectories':	
			# The 'Trajectory' window functions exactly as the field window.
			# So you may therefore display differents things on two windows
			return (('s1',(40 + 20*(((1+self.StepId) % 4)//2), 60 - 20 * ((self.StepId % 4)//2), 'brown', 1, 50, 50, 'brown', 3)),)
			
		#--------------------------------------------------------------------------#
		# Evolving curves	(option C when starting Evolife)                       #
		#--------------------------------------------------------------------------#
		elif Slot == 'PlotOrders':	
			# Should return a list of plot orders for drawing dynamic curves 
			#(typically curves showing the evolution of some values)
			# Each plot order has the form:  (Colour, (x,y))
			return [('blue', (self.StepId, self.StepId + randint(0,20) ))]	# draws a noisy line increasing at each time step
		
		elif Slot == 'CurveNames':	return [('blue', 'time step + noise')]	# provides a legend
		
		#--------------------------------------------------------------------------#
		# Displaying images	(option F when starting Evolife)                       #
		#--------------------------------------------------------------------------#
		elif Slot == 'Image':
			# Should return the path to an image file
			return None

		#--------------------------------------------------------------------------#
		# Displaying patterns (option T when starting Evolife)                     #
		# (same as slot 'Image', but for the 'Trajectories' window)                  #
		#--------------------------------------------------------------------------#
		elif Slot == 'Pattern':
			# Should return the path to an image file
			return None
			
		else:   return EO.Generic_Observer.get_info(self, Slot)	# basic behaviour; 


class Individual:
	"   class Individual: defines what an individual consists of "
	def __init__(self, IdNb):
		self.id = "A%d" % IdNb	# Identity number
		if IdNb % 2:	self.colour = 'red'
		else:			self.colour = 'blue'	# available colours are displayed at the right of Evolife curve display, from bottom up.
		self.Location = (IdNb, randint(0, IdNb), self.colour, 5)
	
	def move(self):	
		# Just a brownian vertical movement
		# self.Location = (self.Location[0], max(0, self.Location[1] + randint(-1, 1)), self.colour, 5)	# 5 == size of blobs
		self.Location = (self.Location[0], max(0, self.Location[1] + randint(-1, 1)), self.colour, 5, 50, 0, 'grey', 1)	# 5 == size of blobs
		
class Population:
	" defines the population of agents "

	def __init__(self, NbAgents, Observer):
		" creates a population of agents "
		self.Pop = [Individual(IdNb) for IdNb in range(NbAgents)]
		self.Obs = Observer
		self.Obs.Positions = self.positions()
		self.Obs.NbAgents = NbAgents
				 
	def positions(self):
		return dict([(A.id, A.Location) for A in self.Pop])
		
	################################################	  
	# This function is run at each simulation step #	
	################################################	  
	def one_year(self):
		self.Obs.StepId += 1	# udates simulation step in Observer
		for agent in self.Pop:
			agent.move()		# whatever behaviour one wants to simulate
		self.Obs.Positions = self.positions()	# let Observer know position changes
		return True




		
def Start():
	Obs = Observer()   # Observer contains statistics
	Obs.setOutputDir('___Results')	# curves, average values and screenshots will be stored there
	Obs.recordInfo('BackGround', 'yellow')	# windows will have this background by default
											# BackGround could be 'yellow' or 'toto.jpg' or '11' or '#F0B554'
	Obs.recordInfo('CurvesWallpaper', '../QtGraphics/EvolifeBG.png')
	Obs.recordInfo('TrajectoriesWallpaper', '../QtGraphics/EvolifeBG.png')
	Obs.recordInfo('DefaultViews',	['Field', 'Trajectories'])	# Evolife should start with that window open
	# Obs.recordInfo('DefaultViews',	['Field'])	# Evolife should start with that window open
	Pop = Population(Gbl.NbAgents, Obs)   # population of agents

	#--------------------------------------------------------------------------#
	# Start: launching Evolife window system                                   #
	#--------------------------------------------------------------------------#
	# Start(Callback function, Observer, Capabilities, Options)
	# Callback function:	this locally defined function is called by Evolife at each time step
	# Observer:	locally defined 
	# Capabilities: string containing any of the following letters
	#	C = Curves 
	#	F = Field (2D seasonal display) (excludes R)
	#	G = Genome display
	#	L = Log Terminal (not implemented)
	#	N = social network display
	#	P = Photo (screenshot)
	#	R = Region (2D ongoing display) (excludes F)
	#	T = Trajectory display (requires additional implementation in each case)

	EW.Start(
		Pop.one_year, 
		Obs, 
		Capabilities='FCPT'
		)
	


if __name__ == "__main__":
		print __doc__
		Start()
		print "Bye......."
		sleep(1.1)	


__author__ = 'Dessalles'
