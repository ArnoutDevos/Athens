##############################################################################
# SocialNetwork   www.dessalles.fr/Evolife              Jean-Louis Dessalles #
#                 Telecom ParisTech  2014                www.dessalles.fr    #
##############################################################################


""" Study of the role of signalling in the emergence of social networks:
Individuals must find a compromise between efforts devoted to signalling
and the probability of attracting followers.
Links are supposed to be symmetrical.
"""


#import math
from time import sleep
from random import randint, random, sample, choice

import sys
import os.path
sys.path.append('../../..')

import Evolife.QtGraphics.Evolife_Window 
import Evolife.Scenarii.Parameters 
from Evolife.Tools.Tools import boost, percent, LimitedMemory
from Evolife.Ecology.Observer import Experiment_Observer
from Evolife.Ecology.Alliances import Alliances
from Evolife.Scenarii.S_Signalling import Interactions


# Global elements
class Global(object):
	def __init__(self):
		# General functions
		self.Closer = lambda x, Target, Attractiveness: ((100.0 - Attractiveness) * x + Attractiveness * Target) / 100
			# pushes x towards Target
		self.Perturbate = lambda x, Amplitude: x + (2 * random() - 1) * Amplitude
			# mutation function
		self.Limitate = lambda x, Min, Max: min(max(x,Min), Max)
			# keeps x within limits
		self.Decrease = lambda x, MaxX, MinY: (100 - x * ((100.0 - MinY)/ MaxX))
			# linear decreasing function between 100 and MinY

		# Parameter values
		self.Parameters = Evolife.Scenarii.Parameters.Parameters('_Params.evo')
		self.ScenarioName = self.Parameters.Parameter('ScenarioName')
		# Definition of interaction
		self.Interactions = None	# to be overloaded

	def Dump_(self, PopDump, ResultFileName, Verbose=False):
		""" Saves parameter values, then agents' investment in signalling, then agents' distance to best friend
		"""
		if Verbose:
			print "Saving data to %s.*" % ResultFileName
		SNResultFile = open(ResultFileName + '.res', 'a')
		SNResultFile.write('\n' + "\t".join(PopDump('SignalInvestment')))	  
		SNResultFile.write('\n' + "\t".join(PopDump('DistanceToBestFriend')))	  
		SNResultFile.close()


Gbl = Global()



if Gbl.Parameters.Parameter('BatchMode'):
	boost()   # A technical trick that sometimes provides impressive speeding up
else:
	print boost()   # A technical trick that sometimes provides impressive speeding up
	print "Interaction scenario: %s\n\n" % Gbl.ScenarioName



class SNetwork_Observer(Experiment_Observer):
	" Stores some global observation and display variables "

	def __init__(self, Parameters):
		Experiment_Observer.__init__(self, Parameters)
		#additional parameters	  
		self.NbAgents = Parameters.Parameter('NbAgents')	 # total number of agents
		self.Positions = dict() # position of agents, for display
		self.Alliances = []
		
	def get_data(self, Slot):
		if Slot == 'Positions':	return self.Positions.items()	
		elif Slot == 'Network': return self.Alliances
		else:   return Experiment_Observer.get_data(self, Slot)

	def get_info(self, Slot):
		if Slot == 'PlotOrders':	return [(3, (self.StepId, self.StepId))]
		elif Slot == 'CurveNames':	return [(3, 'year')]
		else:   return Experiment_Observer.get_info(self, Slot)

if Gbl.ScenarioName== 'Grooming':
	class Inter_Actions(Interactions):
		""" A few functions used to rule interactions between agents
			Inherits from the signalling scenario """
		
		def __init__(self, ParameterSet):
			Interactions.__init__(self, ParameterSet, self.LCompetence)
			self.Transparency = (ParameterSet('Transparency') != 0)
	##		self.SaturateTable = []
	##		self.Saturate(0, Init=True)
	##		self.Parameter = Parameter  # dict of (parameter, value) pairs
	##		self.FCompetence = self.Competence  # Competence function that computes the competence of individuals


	##	def Saturate(self, x, Init=False):
	##		if Init:
	##			for ii in range(100):
	##				self.SaturateTable.append(math.sin((math.pi * ii) / 200))
	##		return self.SaturateTable[int(x)]
		
		def LCompetence(self, Indiv, Apparent=False):
			if self.Transparency:
				return Indiv.Competence(Apparent=False)
			else:
				return Indiv.Competence(Apparent)

	Gbl.InterActions = Inter_Actions(Gbl.Parameters.Parameter)
	
class Individual(Alliances):
	"   class Individual: defines what an individual consists of "

	def __init__(self, IdNb, NbAgents):
		self.id = "A%d" % IdNb	# Identity number
		Alliances.__init__(self,Gbl.Parameters.Parameter('MaxGurus'),Gbl.Parameters.Parameter('MaxFollowers'))
		self.Quality = (100 * IdNb) / NbAgents # quality may be displayed
		self.Renew()
		self.Age = randint(0, Gbl.Parameters.Parameter('AgeMax'))
		self.update()

	def Renew(self):
		self.SignalInvestment = randint(0,100)   # propensity to display quality
		self.Benefits = LimitedMemory(Gbl.Parameters.Parameter('MemorySpan'))  # memory of past benefits
		self.Points = 0
		self.Risk = 0
		self.Age = 0
		self.detach()
		self.update()

	def update(self):
		#if self.adult():
		#Colour = 21 -(10*self.Age)/Gbl.Parameters.Parameter('AgeMax')
		#if not self.best_friend_symmetry(): Colour = 6
		if self.adult():
			Colour = 'brown'
		else:
			Colour = 'pink'
		# self.Signal = percent(self.SignalInvestment * self.Quality)
		#self.Signal = Gbl.InterActions.FCompetence(self, Apparent=True)
		self.Signal = self.Competence(Apparent=True)
		# self.Position = (self.id, self.SignalInvestment, Colour)
		self.Position = (self.Quality, self.Signal+1, Colour, 8)	# 8 == size of blob in display
		if Gbl.Parameters.Parameter('Lines') and self.best_friend():
			self.Position += (self.best_friend().Quality, self.best_friend().Signal, 21, 1)
			
	def Competence(self, Apparent=False):
		" returns the actual quality of an individual or its displayed version "
		BC = Gbl.Parameters.Parameter('BottomCompetence')
		Comp = percent(100-BC) * self.Quality + BC
		# Comp = BC + self.Quality
		VisibleCompetence = percent(Comp * self.SignalInvestment)	   
		if Apparent:
			return Gbl.Limitate(VisibleCompetence, 0, 100)
		else:
			# return Comp
			return self.Quality

	def imitate(self, models, Attractiveness):
		" the individual moves its investment closer to the models' investment "
		TrueModels = [m for m in models if m.adult()]
		if TrueModels == []:
			return
		ModelValues = map(lambda x: x.SignalInvestment, TrueModels)
		Avg = float(sum(ModelValues)) / len(ModelValues)
		self.SignalInvestment = Gbl.Closer(self.SignalInvestment, Avg, Attractiveness)
		self.update()

	def explore(self, Speed, Conservatism):
		"   the individual changes its investment to try to get more points "
		# retrieve the best solution so far
		#print self.Benefits.past
		if self.Benefits.Length() > 0:
			Best = max(self.Benefits.retrieve(), key = lambda x: x[1])[0]
		else:
			Best = self.SignalInvestment
		Target = Gbl.Limitate(Gbl.Perturbate(Best, Speed),0,100)
		self.SignalInvestment = Gbl.Closer(Target, self.SignalInvestment, Conservatism)
		self.update()

	def adult(self):
		return self.Age > percent(Gbl.Parameters.Parameter('AgeMax') * Gbl.Parameters.Parameter('Infancy'))

	def Learns(self, neighbours, Children):
		""" Learns by randomly changing current value.
			Starting point depends on previous success and on neighbours.
			If 'Children' is true, perturbation is larger for children 
		"""
		if self.Age > Gbl.Parameters.Parameter('AgeMax'):
			self.Renew()
		self.Age += 1
		self.imitate(neighbours, Gbl.Parameters.Parameter('ImitationStrength')) # gets SignalInvestment closer to neighbours' values
		if Children and not self.adult():
			# still a kid
			LearningSpeed = Gbl.Decrease(self.Age, Gbl.Parameters.Parameter('Infancy'), Gbl.Parameters.Parameter('LearningSpeed'))
		else:
			LearningSpeed = Gbl.Parameters.Parameter('LearningSpeed')
		self.explore(LearningSpeed, Gbl.Parameters.Parameter('LearningConservatism'))	# compromise between current value and a perturbation of past best value

	def Interact(self, Signallers):
		if Signallers == []:	return
		if Gbl.ScenarioName == 'CostlySignal':
			# The agent chooses the best available Signaller from a sample.
			OldFriend = self.best_friend()
			Signallers.sort(key=lambda S: S.Signal, reverse=True)
			for Signaller in Signallers:
				if Signaller == self:	continue
				if OldFriend and OldFriend.Signal >= Signaller.Signal:
					break	# no available interesting signaller
				if Signaller.followers.accepts(0) >= 0:
					# cool! Self accepted as fan by Signaller.
					if OldFriend is not None and OldFriend != Signaller:
						self.quit_(OldFriend)
					self.new_friend(Signaller, Signaller.Signal)
					break
		elif Gbl.ScenarioName == 'Grooming':
			Gbl.InterActions.groom(self, Signallers[0])

	def assessment(self):
		self.Points -= percent(Gbl.Parameters.Parameter('SignallingCost') * self.SignalInvestment)
		if Gbl.ScenarioName == 'CostlySignal':
			if self.best_friend() is not None:
				self.best_friend().Points += Gbl.Parameters.Parameter('JoiningBonus')
		elif Gbl.ScenarioName == 'Grooming':
			self.restore_symmetry() # Checks that self is its friends' friend
			for (Rank,Friend) in enumerate(self.friends()):
				AgentSocialOffer = Gbl.InterActions.SocialOffer(Gbl.InterActions.FCompetence(self,
																					Apparent=False),
															Rank, self.nbFriends())
				Friend.Risk = percent(Friend.Risk * (100 - percent(Gbl.Parameters.Parameter('CompetenceImpact')
																		   * AgentSocialOffer)))
				
	def wins(self):
		"   stores a benefit	"
		if Gbl.ScenarioName == 'Grooming':
			self.Points += 100 - self.Risk
			# self.Points = Gbl.InterActions.Saturate(self.Points)
		if self.Benefits.last() is not None and self.Benefits.last()[0] == self.SignalInvestment:
			# if SignalInvestment has not changed, replace last benefit
			Previous = self.Benefits.pull()[1]
			self.Benefits.push((self.SignalInvestment, (Previous + self.Points) / 2))
		else:
			self.Benefits.push((self.SignalInvestment, self.Points))
	
	def __repr__(self):
		return "%s[%0.1f]" % (self.id, self.Signal)
		
class Population(object):
	" defines the population of agents "

	def __init__(self, NbAgents, Observer):
		" creates a population of swallow agents "
		self.Pop = [Individual(IdNb, NbAgents) for IdNb in range(NbAgents)]
		self.PopSize = NbAgents
		self.Obs = Observer
		self.Obs.Positions = self.positions()
		self.Obs.NbAgents = self.PopSize
				 
	def positions(self):
		return dict([(A.id, A.Position) for A in self.Pop])

	def neighbours(self, Agent):
		" Returns a list of neighbours for an agent "
		AgentQualityRank = self.Pop.index(Agent)
		return [self.Pop[NBhood] for NBhood in [AgentQualityRank - 1, AgentQualityRank + 1]
				if NBhood >= 0 and NBhood < self.PopSize]
		  
	def One_Run(self):
		self.Obs.StepId += 1
		##################### PROVISOIRE ############################
##		if self.Obs.StepId == 3000:
##			Gbl.Parameters.Params['LearningSpeed'] = 11
		##################### PROVISOIRE ############################
		for Run in range(self.PopSize):
			self.One_Step()
		return True

	def One_Step(self):
		" Interactions take place, then learning "
		# This procedure is repeatedly called by the simulation thread
		# Social interactions
		for agent in self.Pop:
			agent.Points = 0
			agent.Risk = 100	# maximum risk
			# agent.lessening_friendship()	# eroding past gurus performances
##		for agent in self.Pop:
##			# first interact with current friends
##			for Partner in agent.friends():
##				agent.Interact([Partner])
		for Run in range(Gbl.Parameters.Parameter('NbRuns')):
			if Gbl.ScenarioName == 'CostlySignal':
				Fan = choice(self.Pop)
				# Fan chooses from a sample
				Fan.Interact(sample(self.Pop, Gbl.Parameters.Parameter('SampleSize')))
			elif Gbl.ScenarioName == 'Grooming':
				(Player, Partner) = sample(self.Pop, 2)
				Partner.Interact([Player])
		for agent in self.Pop:
			# one looks whether interactions were profitable
			agent.assessment()
		for agent in self.Pop:
			agent.wins()	# memorizes success
			#if self.Obs.Visible() and agent.Quality == 99:
			#	print agent, '>>>', agent.Benefits, '{%d}' % agent.Points
			#	print "%s --> %s" % (agent, agent.Benefits)
		# Learning occurs for one agent every NbRuns interactions
		Learner = choice(self.Pop)
		Children = self.Obs.StepId < percent(self.Obs.TimeLimit * Gbl.Parameters.Parameter('LearnHorizon'))
		Learner.Learns(self.neighbours(Learner), Children)
		# Learner should re-interact with its followers
		for Partner in Learner.followers.names():
			Partner.Interact([Learner])
		self.Obs.Positions[Learner.id] = Learner.Position
		if self.Obs.Visible():
##			self.Obs.Alliances = [(agent.id,[(agent.best_friend().id,0)])
##									  for agent in self.Pop if agent.best_friend() is not None]
##			self.Obs.Alliances += [(agent.id,[])
##									  for agent in self.Pop if agent.best_friend() == None]
			self.Obs.Alliances = [(agent.id, [T.id for T in Alliances.signature(agent)]) for agent in self.Pop]
		return True	 # This value is forwarded to "ReturnFromThread"

	def Dump(self, Slot):
		""" Saving investment in signalling for each adult agent
			and then distance to best friend for each adult agent having a best friend
		"""
		if Slot == 'SignalInvestment':
			#D = [(agent.Quality, "%2.03f" % agent.SignalInvestment) for agent in self.Pop if agent.adult()]
			#D += [(agent.Quality, " ") for agent in self.Pop if not agent.adult()]
			D = [(agent.Quality, "%2.03f" % agent.SignalInvestment) for agent in self.Pop]
		if Slot == 'DistanceToBestFriend':
			D = [(agent.Quality, "%d" % abs(agent.best_friend().Quality - agent.Quality)) for agent in self.Pop if agent.adult() and agent.best_friend() is not None]
			D += [(agent.Quality, " ") for agent in self.Pop if agent.best_friend() == None or not agent.adult()]
		return [Slot] + [d[1] for d in sorted(D, key=lambda x: x[0])]

		


		
def Start(BatchMode=False):
	Observer = SNetwork_Observer(Gbl.Parameters)   # Observer contains statistics
	Observer.setOutputDir('___Results')
	Observer.recordInfo('DefaultViews',	['Field', 'Network'])	# Evolife should start with that window open
	Pop = Population(Observer.NbAgents, Observer)   # population of agents
	if BatchMode :
		for Step in range(Gbl.Parameters.Parameter('TimeLimit')):
			#print '.',
			Pop.One_Run()
			if os.path.exists('stop'):
				break
		# writing header to result file
		open(Observer.get_info('ResultFile')+'.res','w').write(Observer.get_info('ResultHeader'))
		Gbl.Dump_(Pop.Dump, Observer.get_info('ResultFile'))
		return

	# Interactive mode
	if Gbl.ScenarioName == 'Grooming':
		print 'RankEffects:', Gbl.InterActions.RankEffects	
	" launching window "
	Evolife.QtGraphics.Evolife_Window.Start(Pop.One_Run, Observer, Capabilities='FNCP', Options=[('BackGround','white')])

	Gbl.Dump_(Pop.Dump, Observer.get_info('ResultFile'), Verbose=True)
	


if __name__ == "__main__":

	BatchMode = Gbl.Parameters.Parameter('BatchMode')

	if BatchMode:
		Start(BatchMode=True)
	else:
		print __doc__
		Start(BatchMode=False)
		print "Bye......."

		sleep(2.1)	


__author__ = 'Dessalles'
