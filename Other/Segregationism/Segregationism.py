#!/usr/bin/env python
##############################################################################
# EVOLIFE  www.dessalles.fr/Evolife                    Jean-Louis Dessalles  #
#            Telecom ParisTech  2014                       www.dessalles.fr  #
##############################################################################

##############################################################################
# Segregationism (after Thomas Schelling's work)                             #
##############################################################################

#  Basic Version

""" Emergence of segregationism:
Though individual agents show only slight preference for being surrounded by similar agent, homogeneous patches emerge.
"""
		# Thomas Schelling (1971) studied the dynamics of residential segregation
		# to elucidate the conditions under which individual decisions about where to live
		# will interact to produce neighbourhoods that are segregated by race.
		# His model shows that this can occur even though individuals do not act
		# in a coordinated fashion to bring about these segregated outcomes.
		# Schelling proposed a prototype model in which individual agents are of two types,
		# say red and blue, and are placed randomly on the squares of a checkerboard.
		# The neighbourhood of an agent is defined to be the eight squares adjoining his location.
		# Each agent has preferences over the composition of his neighbourhood,
		# defined as the proportion of reds and blues. In each period, the most dissatisfied
		# agent moves to an empty square provided a square is available that he prefers
		# to his current location. The process continues until no one wants to move.

		# The typical outcome is a highly segregated state, although nobody actually
		# prefers segregation to integration. 


import sys
sys.path.append('..')
sys.path.append('../../..')

import Evolife.Ecology.Observer				as EO
import Evolife.Scenarii.Parameters 			as EPar
import Evolife.Scenarii.Default_Scenario	as ED
import Evolife.QtGraphics.Evolife_Window 	as EW
import Evolife.QtGraphics.Curves			as EC
import Evolife.Ecology.Individual			as EI
import Evolife.Ecology.Group				as EG
import Evolife.Ecology.Population			as EP
import Landscapes

from Evolife.Tools.Tools import boost
print boost()   # A technical trick that sometimes provides impressive speeding up

	
import random


class Scenario(EPar.Parameters):
	def __init__(self):
		# Parameter values
		EPar.Parameters.__init__(self, CfgFile='_Params.evo')	# loads parameters from configuration file
		#############################
		# Global variables		    #
		#############################
		AvailableColours = ['red', 'blue', 'brown', 'yellow', 7] + range(8, 21)
		self.Colours = AvailableColours[:self.Parameter('NbColours')]	# corresponds to Evolife colours
		self.addParameter('NumberOfGroups', self.Parameter('NbColours'))	# useful to create coloured groups

		
class Settling_Observer(EO.Observer):
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
		if Slot == 'Positions':
			CC = self.CurrentChanges
			self.CurrentChanges = []
			return tuple(CC)
		else:	return EO.Observer.get_data(self, Slot)

class Individual(EI.Individual):
	""" Defines individual agents
	"""
	def __init__(self, Scenario, ID=None, Newborn=True):
		EI.Individual.__init__(self, Scenario, ID=ID, Newborn=Newborn)
		self.Colour = Gbl.Colours[0]	# just to initialize
		# print 'creating', self.ID
		self.moves()	# gets a location

	def setColour(self, Colour):	
		self.erase()	# moves away from Land
		self.Colour = Colour
		self.locate(self.location, Erase=False)

	def locate(self, NewPosition, Erase=True):
		" place individual at a specific location on the ground "
		# print 'locating', self, 'at', NewPosition
		if NewPosition and not Land.Modify(NewPosition, self.Colour, check=True): # new position on Land
			return False		 # NewPosition is not available  
		if Erase and self.location and not Land.Modify(self.location, None):
			# erasing previous position on Land
			print 'Error, agent %s badly placed' % self.ID
		self.location = NewPosition
		Observer.record((self.ID, self.location + (self.Colour, self.Scenario.Parameter('DotSize')))) # for ongoing display
		return True

	def erase(self):
		" erase individual from the ground "
		# print 'erasing', self
		if self.location:
			if not Land.Modify(self.location, None):	# erase on Land
				print 'Error, agent %s was badly placed' % self.ID, self.location
			# sending negative colour to display to erase the agent
			NColour = '-' + str(self.Colour)
			Observer.record((self.ID, self.location + (NColour, self.Scenario.Parameter('DotSize'))))

	def decisionToMove(self):
		if self.location is None:	return False # may happen if there is no room left	
		Statistics = Land.InspectNeighbourhood(self.location, self.Scenario.Parameter('NeighbourhoodRadius'))	# Dictionary of colours
		Same = Statistics[self.Colour]
		Different = sum([Statistics[C] for C in self.Scenario.Colours if C != self.Colour])	

						############# (TO BE MODIFIED) #####################################
						############# (TO BE MODIFIED) #####################################
						# use:  self.Scenario.Parameter('Tolerance')
		return True
		return False
	
	def moves(self, Position=None):
		# print 'moving', self
		if Position:
			return self.locate(Position)
		else:
			# pick a random location and go there (TO BE MODIFIED)
			for ii in range(10): # should work at first attempt most of the time
				Landing = Land.randomPosition(Content=None, check=True)	# selects an empty cell
				if Landing and self.locate(Landing):
					return True
				elif ii == 0:	Land.statistics()   # need to update list of available positions
			print "Unable to move to", Position,
			return False

	def dies(self):
		" get off from the Land when dying "
		self.erase()
		EI.Individual.dies(self)
			
	def __repr__(self):
		return "(%s,%s) --> " % (self.ID, self.Colour) + str(self.location)

class Group(EG.Group):
	# The group is a container for individuals.
	# Individuals are stored in self.members

	def __init__(self, Scenario, ID=1, Size=100):
		EG.Group.__init__(self, Scenario, ID, Size)
		
	def setColour(self, Colour):
		for member in self.members:	member.setColour(Colour)	# gives colour to all members
		
	def createIndividual(self, ID=None, Newborn=True):
		# calling local class 'Individual'
		Indiv = Individual(self.Scenario, ID=self.free_ID(), Newborn=Newborn)
		# Individual creation may fail if there is no room left
		if Indiv.location == None:	return None
		return Indiv
		
			
	
class Population(EP.Population):
	" defines the population of agents "
	
	def __init__(self, Scenario, Observer):
		" creates a population of agents "
		EP.Population.__init__(self, Scenario, Observer)
		self.Colours = self.Scenario.Colours
		print self.Colours
		# print "Existing colours: %s" % self.Colours
		for Colour in self.Colours:
			print "creating %s agents" % Colour
			# ColourNbr = EC.EvolifeColourID(Colour)[0]	# colour number according to Evolife
			# individuals are created with the colour given as ID of their group
			self.groups[self.Colours.index(Colour)].setColour(Colour)
		print "population size: %d" % self.popSize
		self.Moves = 0  # counts the number of times agents have moved
		self.CallsSinceLastMove = 0  # counts the number of times agents were proposed to move since last actual move

	def createGroup(self, ID=0, Size=0):
		return Group(self.Scenario, ID=ID, Size=Size)

	def One_Decision(self):
		""" This function is repeatedly called by the simulation thread.
			One agent is randomly chosen and decides what it does
		"""

		EP.Population.one_year(self)	# performs statistics

		agent = self.selectIndividual()	# agent who will play the game	
		# print agent.ID, 'about to move'
		self.CallsSinceLastMove += 1
		if agent.decisionToMove() and agent.moves():
			self.Moves += 1
			self.CallsSinceLastMove = 0
		if self.CallsSinceLastMove > 10 * self.popSize:
			return False	# situation is probably stable
		return True
			   

if __name__ == "__main__":
	print __doc__

	
	#############################
	# Global objects			#
	#############################
	Gbl = Scenario()
	Observer = Settling_Observer(Gbl)	  # Observer contains statistics
	Land = Landscapes.Landscape(Gbl.Parameter('LandSize'))	  # logical settlement grid
	Land.setAdmissible(Gbl.Colours)
	Pop = Population(Gbl, Observer)   
	
	# Observer.recordInfo('BackGround', 'white')
	Observer.recordInfo('FieldWallpaper', 'white')
	
	EW.Start(Pop.One_Decision, Observer, Capabilities='RP')

	print "Bye......."
	
__author__ = 'Dessalles'
