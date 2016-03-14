#!/usr/bin/env python
##############################################################################
# SocialSimulation   www.dessalles.fr/Evolife           Jean-Louis Dessalles #
#                 Telecom ParisTech  2014                www.dessalles.fr    #
##############################################################################


""" A basic framework to run social simulations
"""


from time import sleep
from random import sample, randint

import sys
import os.path
sys.path.append('../../..')	# to include path to Evolife

import Evolife.QtGraphics.Evolife_Window 
import Evolife.QtGraphics.Evolife_Batch 
import Evolife.Scenarii.Parameters 
import Evolife.Tools.Tools as ET
import Evolife.Ecology.Observer as EO
import Evolife.Ecology.Alliances as EA
import Evolife.Ecology.Learner as EL



class Social_Observer(EO.Experiment_Observer):
	" Stores some global observation and display variables "

	def __init__(self, Parameters):
		EO.Experiment_Observer.__init__(self, Parameters)
		#additional parameters	  
		self.Positions = dict() # position of agents, for display
		self.Alliances = []		# social links, for display
		self.curve('year', 'blue', 'year')	# default curve: just draw time
		
	def get_data(self, Slot):
		if Slot == 'Positions':	return self.Positions.items()	# returns stored positions
		elif Slot == 'Network':	return self.Alliances			# returns stored links
		else:   return EO.Experiment_Observer.get_data(self, Slot)

	def get_info(self, Slot):
		if Slot == 'PlotOrders':	return [(self.curve('year'), (self.StepId, self.StepId))]
		elif Slot == 'ResultOffset': return self.Parameter('DumpStart') # minimum x-value when computing curve average 
		else:   return EO.Experiment_Observer.get_info(self, Slot)

	def hot_phase(self):
		return self.StepId < self.TimeLimit * self.Parameter('LearnHorizon') / 100.0

	def season(self):	EO.Experiment_Observer.season(self)		# merely increments self.StepId
		
class Social_Individual(EA.Follower, EL.Learner):
	"   A social individual has friends and can learn "

	def __init__(self, IdNb, features={}, maxQuality=100, parameters=None):
		if parameters: 	self.Param = parameters.Param
		else:	self.Param = None	# but this will provoke an error
		self.id = "A%d" % IdNb	# Identity number
		if self.Param('SocialSymmetry'):
			EA.Follower.__init__(self, self.Param('MaxFriends'))
		else:
			EA.Follower.__init__(self, self.Param('MaxFriends'), self.Param('MaxFollowers'))
		self.Quality = (100.0 * IdNb) / maxQuality # quality may be displayed
		# Learnable features
		EL.Learner.__init__(self, features, MemorySpan=self.Param('MemorySpan'), AgeMax=self.Param('AgeMax'), 
							Infancy=self.Param('Infancy'), Imitation=self.Param('ImitationStrength'), 
							Speed=self.Param('LearningSpeed'), Conservatism=self.Param('LearningConservatism'), toric=self.Param('Toric'))
		self.Points = 0	# stores current performance
		self.update()

	def Reset(self):	# called by Learner when born again
		self.forgetAll()	# erase friendships
		EL.Learner.Reset(self)
		
	def update(self, infancy=True):
		"	updates values for display "
		Colour = 'green%d' % int(1 + 10 * (1 - float(self.Age)/(1+self.Param('AgeMax'))))
		if infancy and not self.adult():	Colour = 'red'
		if self.Features:	y = self.Features[self.Features.keys()[0]]
		else:	y = 17
		self.Position = (self.Quality, Colour)


	def Interact(self, Partner):	
		pass	# to be overloaded
		return True

	def assessment(self):
		" Social benefit from having friends "
		pass		# to be overloaded
		return self.Points	
		
	def __str__(self):
		return "%s[%s]" % (self.id, str(self.Features))
		
class Social_Population:
	" defines a population of interacting agents "

	def __init__(self, parameters, NbAgents, Observer, IndividualClass=None):
		" creates a population of agents "
		if IndividualClass is None:	IndividualClass = Social_Individual
		self.Pop = [IndividualClass(IdNb, maxQuality=NbAgents) for IdNb in range(NbAgents)]
		self.PopSize = NbAgents
		self.Obs = Observer
		self.Obs.Positions = self.positions()
		self.Param = parameters.Param
				 
	def positions(self):
		return dict([(A.id, A.Position) for A in self.Pop])

	def neighbours(self, Agent):
		" Returns a list of neighbours for an agent "
		AgentQualityRank = self.Pop.index(Agent)
		return [self.Pop[NBhood] for NBhood in [AgentQualityRank - 1, AgentQualityRank + 1]
				if NBhood >= 0 and NBhood < self.PopSize]
		  
	def display(self, averages=False):
		if self.Obs.Visible():	# Statistics for display
			for agent in self.Pop:
				agent.update(infancy=self.Obs.hot_phase())	# update position for display
				self.Obs.Positions[agent.id] = agent.Position	# Observer stores agent position 
			# Observer stores social links
			self.Obs.Alliances = [(agent.id, [T.id for T in agent.social_signature()]) for agent in self.Pop]
			if averages and self.Pop:
				for F in self.Pop[0].Features:	# First agent's features (considered typical)
					self.Obs.AvgFeatures[F] = float(sum([agent.feature(F) for agent in self.Pop]))/len(self.Pop)

	def season_initialization(self):
		for agent in self.Pop:
			# agent.lessening_friendship()	# eroding past gurus performances
			if self.Param('EraseNetwork'):	agent.forgetAll()
			agent.Points = 100
	
	def interactions(self):
		successful = 0
		for Run in range(self.Param('NbInteractions')):
			(Player, Partner) = sample(self.Pop, 2)
			if Partner.Interact(Player):	successful += 1
		# if randint(1, 100) < 10:
			# print '%d%%' % int((100.0 * successful) / self.Param('NbInteractions')),
	
	def learning(self):
		for agent in self.Pop:	agent.assessment()	# storing current scores (with possible cross-benefits)
		# some agents learn
		Learners = sample(self.Pop, ET.chances(self.Param('LearningProbability')/100.0, len(self.Pop)))	
		for agent in self.Pop:
			agent.wins(agent.Points)	# Stores points for learning
			if agent in Learners:
				if not agent.Learns(self.neighbours(agent), hot=self.Obs.hot_phase()):
					# agent.update()	# this is a newborn - update position for display
					pass
				agent.update()	# update position for display
				# self.Obs.Positions[agent.id] = agent.Position
				
	def One_Run(self, averages=False):
		# This procedure is repeatedly called by the simulation thread
		# ====================
		# Display
		# ====================
		self.Obs.season()	# increments year
		self.display(averages=averages)
		# ====================
		# Interactions
		# ====================
		for Run in range(self.Param('NbRunPerYear')):	
			self.season_initialization()
			self.interactions()
			self.learning()
		return True	# This value is forwarded to "ReturnFromThread"

		
		
def Start(Population, Observer, BatchMode=False):
	
	####################
	# Batch mode
	####################
	if BatchMode :
		# # # # for Step in range(Gbl.Param('TimeLimit')):
			# # # # #print '.',
			# # # # Pop.One_Run()
			# # # # if os.path.exists('stop'):	break
		Evolife.QtGraphics.Evolife_Batch.Start(Population.One_Run, Observer)
		# writing header to result file
		open(Observer.get_info('ResultFile')+'.res','w').write(Observer.get_info('ResultHeader'))
		return

	####################
	# Interactive mode
	####################
	print __doc__
	" launching window "
	# Evolife.QtGraphics.Evolife_Window.Start(Pop.One_Run, Observer, Capabilities='FNCP', Options=[('BackGround','grey')])
	Evolife.QtGraphics.Evolife_Window.Start(Population.One_Run, Observer, Capabilities='FNCP', Options=[('BackGround','lightblue')])

	print "Bye......."
	sleep(2.1)	
	return
	


if __name__ == "__main__":
	ET.boost()   # A technical trick that sometimes provides impressive speeding up with Python up to 2.6
	Observer = Social_Observer(Gbl.Parameters)   # Observer contains statistics
	Observer.setOutputDir('___Results')
	Observer.recordInfo('DefaultViews',	['Field', 'Network'])	# Evolife should start with that window open
	Pop = Population(Gbl.Param('NbAgents'), Observer)   # population of agents
	Start(BatchMode=Gbl.Param('BatchMode'))



__author__ = 'Dessalles'
