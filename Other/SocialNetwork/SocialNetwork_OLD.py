#!/usr/bin/env python
##############################################################################
# SocialNetwork                                         Jean-Louis Dessalles #
#              Telecom ParisTech  2014                   www.dessalles.fr    #
##############################################################################


""" Study of the role of signalling in the emergence of social networks:
Individuals must find a compromise between efforts devoted to signalling
and the probability of attracting followers.
Links are supposed to be symmetrical.
"""

try:
	import psyco	# A technical trick - look at http://psyco.sourceforge.net/
					#(somewhat magical, but sometimes provides impressive speeding up)
	psyco.profile()
	#print '(boosting Python)\n'
except:
	#print '(slow Python)\n'
	pass

#import math
from time import sleep, strftime
from random import randint, random, sample, choice

if __name__ == "__main__":
	# Putting Evolife into the path (supposing we are in the same directory tree)
	import sys, os, os.path
	for R in os.walk(os.path.abspath('.')[0:os.path.abspath('.').find('Evo')]):
		if os.path.exists(os.path.join(R[0],'Evolife','__init__.py')):
			sys.path.append(os.path.abspath(R[0]))
			break

##global MyScenario
##import Evolife.Scenarii.MyScenario
##MyScenario = Evolife.Scenarii.MyScenario.InstantiateScenario(Configuration, '', '')


from Evolife.Scenarii.Parameters import Parameters
MyScenario = Parameters('SNetwork_.evo')

from Evolife.Ecology.Observer import Generic_Observer
from Evolife.Ecology.Alliances import Alliances
from Evolife.Scenarii.S_Signalling import Interactions
from Evolife.Tools.Tools import percent, LimitedMemory

if MyScenario.Parameter('BatchMode') == 0:
	import Tkinter
	from Evolife.Graphics.Plot_Area import Ground
	from Evolife.Graphics.Evolife_Graphic import Alliance_Area
	from Evolife.Graphics.Evolife_Window import Satellite_window, Generic_Main_Frame


# General functions
Closer = lambda x, Target, Attractiveness: ((100.0 - Attractiveness) * x + Attractiveness * Target) / 100
	# pushes x towards Target
Perturbate = lambda x, Amplitude: x + (2 * random() - 1) * Amplitude
	# mutation function
Limitate = lambda x, Min, Max: min(max(x,Min), Max)
	# keeps x within limits
Decrease = lambda x, MaxX, MinY: (100 - x * ((100.0 - MinY)/ MaxX))
	# linear decreasing function between 0 and MaxX



class SNetwork_Observer(Generic_Observer):
	" Stores some global observation and display variables "

	def __init__(self, ExperienceName):
		Generic_Observer.__init__(self, 'SNetwork_Observer')

		# overloaded parameters
		self.TimeLimit = MyScenario.Parameter('TimeLimit')
		self.ScenarioName = 'SNetwork'
		self.ExperienceID = ExperienceName
		self.DispPeriod = MyScenario.Parameter('DispPeriod')
		
		#additional parameters	  
		self.NbAgents = MyScenario.Parameter('NbAgents')	 # total number of agents
		self.Positions = [] # position of agents, for display
		self.Alliances = []
		
class Inter_Actions(Interactions):
	" A few functions used to rule interactions between agents "
	
	def __init__(self, Parameter):
		Interactions.__init__(self, Parameter, self.Competence)
##        self.SaturateTable = []
##        self.Saturate(0, Init=True)
##        self.Parameter = Parameter  # dict of (parameter, value) pairs
##        self.FCompetence = self.Competence  # Competence function that computes the competence of individuals

	def Competence(self, Indiv, Apparent=False):
		" returns the actual quality of an individual or its displayed version "
		BC = self.Parameter('BottomCompetence')
		Comp = percent(100-BC) * Indiv.Quality + BC
		# Comp = BC + Indiv.Quality
		VisibleCompetence = percent(Comp * Indiv.SignalInvestment)	   
		if Apparent:
			return Limitate(VisibleCompetence, 0, 100)
		else:
			# return Comp
			return Indiv.Quality

##    def Saturate(self, x, Init=False):
##        if Init:
##            for ii in range(100):
##                self.SaturateTable.append(math.sin((math.pi * ii) / 200))
##        return self.SaturateTable[int(x)]
		

InterActions = Inter_Actions(MyScenario.Parameter)
	
class Individual(Alliances):
	"   class Individual: defines what an individual consists of "

	def __init__(self, IdNb):
		self.id = IdNb	# Identity number
		Alliances.__init__(self,MyScenario.Parameter('MaxGurus'),MyScenario.Parameter('MaxFollowers'))
		self.Quality = IdNb # quality may be displayed
		self.Renew()
		self.Age = randint(0, MyScenario.Parameter('AgeMax'))
		self.update()

	def Renew(self):
		self.SignalInvestment = randint(0,100)   # propensity to display quality
		self.Benefits = LimitedMemory(MyScenario.Parameter('MemorySpan'))  # memory of past benefits
		self.Points = 0
		self.Risk = 0
		self.Age = 0
		self.detach()
		self.update()

	def update(self):
		#if self.adult():
		#Colour = 21 -(10*self.Age)/MyScenario.Parameter('AgeMax')
		#if not self.best_friend_symmetry(): Colour = 6
		if self.adult():
			Colour = 6
		else:
			Colour = 4
		# self.Signal = percent(self.SignalInvestment * self.Quality)
		self.Signal = InterActions.FCompetence(self, Apparent=True)
		# self.Position = (self.id, self.SignalInvestment, Colour)
		self.Position = (self.id, self.Signal, Colour)
		
	def imitate(self, models, Attractiveness):
		" the individual moves its investment closer to the models' investment "
		TrueModels = [m for m in models if m.adult()]
		if TrueModels == []:
			return
		ModelValues = map(lambda x: x.SignalInvestment, TrueModels)
		Avg = float(sum(ModelValues)) / len(ModelValues)
		self.SignalInvestment = Closer(self.SignalInvestment, Avg, Attractiveness)
		self.update()

	def explore(self, Speed, Conservatism):
		"   the individual changes its investment to try to get more points "
		# retrieve the best solution so far
		#print self.Benefits.past
		if self.Benefits.Length() > 0:
			Best = max(self.Benefits.retrieve(), key = lambda x: x[1])[0]
		else:
			Best = self.SignalInvestment
		Target = Limitate(Perturbate(Best, Speed),0,100)
		self.SignalInvestment = Closer(Target, self.SignalInvestment, Conservatism)
		self.update()

	def adult(self):
		return self.Age > percent(MyScenario.Parameter('AgeMax') * MyScenario.Parameter('Infancy'))

	def Learns(self, neighbours, Children):
		""" Learns by randomly changing current value.
			Starting point depends on previous success and on neighbours.
			If 'Children' is true, perturbation is larger for children 
		"""
		if self.Age > MyScenario.Parameter('AgeMax'):
			self.Renew()
		self.Age += 1
		self.imitate(neighbours, MyScenario.Parameter('ImitationStrength'))
		if Children and not self.adult():
			# still a kid
			LearningSpeed = Decrease(self.Age, MyScenario.Parameter('Infancy'), MyScenario.Parameter('LearningSpeed'))
		else:
			LearningSpeed = MyScenario.Parameter('LearningSpeed')
		self.explore(LearningSpeed, MyScenario.Parameter('LearningConservatism'))

	def Interact(self, Signaller):
		# The agent may choose the Signaller as new friend.
		# But she must first check her old friends
		if MyScenario.Parameter('InteractionScenario') == MyScenario.Parameter('CostlySignal'):
			OldFriend = self.best_friend()
			if OldFriend == None \
			   or OldFriend.Signal < Signaller.Signal \
			   or OldFriend == Signaller:   # as Signaller may have changed signal
				if OldFriend != None:
					self.quit_(OldFriend)
				self.new_friend(Signaller, Signaller.Signal)
		elif MyScenario.Parameter('InteractionScenario') == MyScenario.Parameter('Grooming'):
			InterActions.groom(self, Signaller)

	def assessment(self):
		self.Points -= percent(MyScenario.Parameter('SignallingCost') * self.SignalInvestment)
		if MyScenario.Parameter('InteractionScenario') == MyScenario.Parameter('CostlySignal'):
			if self.best_friend() is not None:
				self.best_friend().Points += MyScenario.Parameter('JoiningBonus')
		elif MyScenario.Parameter('InteractionScenario') == MyScenario.Parameter('Grooming'):
			self.restore_symmetry() # Checks that self is its friends' friend
			for (Rank,Friend) in enumerate(self.friends()):
				AgentSocialOffer = InterActions.SocialOffer(InterActions.FCompetence(self,
																					Apparent=False),
															Rank, self.nbFriends())
				Friend.Risk = percent(Friend.Risk * (100 - percent(MyScenario.Parameter('CompetenceImpact')
																		   * AgentSocialOffer)))
				
	def wins(self):
		"   stores a benefit	"
		if MyScenario.Parameter('InteractionScenario') == MyScenario.Parameter('Grooming'):
			self.Points += 100 - self.Risk
			# self.Points = InterActions.Saturate(self.Points)
		if self.Benefits.last() is not None and self.Benefits.last()[0] == self.SignalInvestment:
			# if SignalInvestment has not changed, replace last benefit
			Previous = self.Benefits.pull()[1]
			self.Benefits.push((self.SignalInvestment, (Previous + self.Points) / 2))
		else:
			self.Benefits.push((self.SignalInvestment, self.Points))
	
	def __repr__(self):
		return "I%d[%0.1f]" % (self.id, self.Signal)
		
class Population(object):
	" defines the population of agents "

	def __init__(self, NbAgents, Observer):
		" creates a population of swallow agents "
		self.Pop = [Individual(IdNb) for IdNb in range(NbAgents)]
		self.PopSize = NbAgents
		self.Obs = Observer
		self.Obs.Positions = self.positions()
		self.Obs.NbAgents = self.PopSize
				 
	def positions(self):
		return [A.Position for A in self.Pop]

	def neighbours(self, Agent):
		" Returns a list of neighbours for an agent "
		return [self.Pop[NBhood] for NBhood in [Agent.id - 1, Agent.id + 1]
				if NBhood >= 0 and NBhood < self.PopSize]
		  
	def One_Run(self):
		self.Obs.StepId += 1
		##################### PROVISOIRE ############################
##        if self.Obs.StepId == 3000:
##            MyScenario.Params['LearningSpeed'] = 11
		##################### PROVISOIRE ############################
		for Run in range(self.PopSize):
			self.One_Step()

	def One_Step(self):
		" Interactions take place, then learning "
		# This procedure is repeatedly called by the simulation thread
		# Social interactions
		for agent in self.Pop:
			agent.Points = 0
			agent.Risk = 100	# maximum risk
##        for agent in self.Pop:
##            # first interact with current friends
##            for Partner in agent.friends():
##                agent.Interact(Partner)
		for Run in range(MyScenario.Parameter('NbRuns')):
			(Player, Partner) = sample(self.Pop, 2)
			Partner.Interact(Player)
		for agent in self.Pop:
			# one looks whether interactions were profitable
			agent.assessment()
		for agent in self.Pop:
			agent.wins()	# memorizes success
			#if self.Obs.Visible() and agent.id == 99:
			#	print agent, '>>>', agent.Benefits, '{%d}' % agent.Points
			#	print "%s --> %s" % (agent, agent.Benefits)
		# Learning occurs for one agent every NbRuns interactions
		Learner = choice(self.Pop)
		Children = self.Obs.StepId < percent(self.Obs.TimeLimit * MyScenario.Parameter('LearnHorizon'))
		Learner.Learns(self.neighbours(Learner), Children)
		# Learner should re-interact with its followers
		for Partner in Learner.followers.names():
			Partner.Interact(Learner)
		self.Obs.Positions[Learner.id] = Learner.Position
		if self.Obs.Visible():
##            self.Obs.Alliances = [(agent.id,[(agent.best_friend().id,0)])
##                                      for agent in self.Pop if agent.best_friend() is not None]
##            self.Obs.Alliances += [(agent.id,[])
##                                      for agent in self.Pop if agent.best_friend() == None]
			self.Obs.Alliances = [(agent.id,Alliances.signature(agent)) for agent in self.Pop]
		return 0	 # This value is forwarded to "ReturnFromThread"

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

if MyScenario.Parameter('BatchMode') == 0:
	class SNetwork_Frame(Generic_Main_Frame):

		def __init__(self, Parent, SimulationStep, Observer):
			""" creation of a generic main frame, including control buttons
				and a graphic area for curve display """
			Generic_Main_Frame.__init__(self, Parent, SimulationStep, Observer, Wtitle='Social Network')
			self.Obs = Observer
			self.StartGround()
			self.StartAllianceArea()
			self.LastTimeStep = 0
			self.FilmMode = 0
			
		def StartGround(self):
			""" the ground is used to plot the investment of agents in communication
			"""
			if 'SNetwork' not in self.SWindows:
				# the ground has to be created
				self.CreateSatelliteWindowContainer('SNetwork')
				self.Ground = Ground(self.SWindows['SNetwork'], Toric=False)
				self.Ground.pack(expand=Tkinter.YES,fill=Tkinter.BOTH)
				for AgentIdNb in range(len(self.Obs.Positions)):
					self.Ground.create_agent('S'+str(AgentIdNb),
											 self.Obs.Positions[AgentIdNb])
			#self.SWindows['SNetwork'].focus_force()
			self.SWindows['SNetwork'].Raise()

		def StartAllianceArea(self):
			""" Social links are displayed as lines between two horizontal lines
				where individuals are located (twice)
			"""
			if 'SLinks' not in self.SWindows:
				# the alliance area has to be created
				self.CreateSatelliteWindowContainer('SLinks')
				self.SLinks = Alliance_Area(self.SWindows['SLinks'])
				self.SLinks.pack(expand=Tkinter.YES,fill=Tkinter.BOTH)
				self.Alliance_display()
			# self.SWindows['SLinks'].focus_force()
			self.SWindows['SLinks'].Raise()

		def Alliance_display(self):
			if 'SLinks' in self.SWindows:
				Layout = map(lambda x: (x[0],x), self.Obs.Positions) # Agents' horizontal positions
				if self.Obs.Alliances != []:
					self.SLinks.Alliance_display(Layout, self.Obs.Alliances)
		
		def GroundUpdate(self):
			" creates dots representing agents "
			if 'SNetwork' in self.SWindows:
				for AgentIdNb in range(len(self.Obs.Positions)):
					self.Ground.move_agent('S'+str(AgentIdNb), self.Obs.Positions[AgentIdNb])
			
		def ReturnFromThread(self, AgentIdNb):
			" Actions to be performed when 'SimulationStep' has been performed "
			# this function is called back after each simulation step
			if self.Busy():
				print 'Busy...',
				sleep(3)
				self.Busy(0, Force=True)
				print '<>'
			if self.Obs.Visible():
				self.Busy(1)  # diving deeper into busy mode
				try:
					self.GroundUpdate()
					self.Alliance_display()
				except ValueError:
					print 'Display jam...'
				if self.LastTimeStep != self.Obs.StepId:
					self.plot_area.plot(3, (self.Obs.StepId, self.Obs.StepId))
					self.LastTimeStep = self.Obs.StepId
				sleep(0.5)	# necessary to avoid graphic jam
				if self.FilmMode:
					self.photo()
				self.after_idle(self.Busy,-1)  # one step back from busy mode
			if self.Obs.Over():
				self.End()
				return -1	# Stops the simulation thread
			return 0 # the simulation resumes

		def photo(self):
			self.Ground.photo('SN_communicationInvestment')
			self.SLinks.photo('SN_Links')
			#self.SLinks.Alliance_display(Layout, self.Obs.Alliances, Photo=True)


		def KeyPressed(self, event):
			Generic_Main_Frame.KeyPressed(self, event)
			# Additional key actions 
			if event.char in ['g','f']:
				self.StartGround()
			elif event.char in ['a','A']:
				self.StartAllianceArea()
			elif event.char == 'p':
				self.photo()
			elif event.char == 'P':
				self.FilmMode = 1 - self.FilmMode

		def End(self):
			print 'Dumping...'
			self.photo()
		

def Dump(Pop, ExpeName, LastStep):
	""" Saves parameter values, then agents' investment in signalling, then agents' distance to best friend
	"""
	SNResultFile = open('SN_%s.res' % ExpeName, 'w')
	SNResultFile.write('Date\t' + '\t'.join(MyScenario.RelevantParamNames()) + '\tLastStep')
	#SNResultFile.write('\tLastStep')
	SNResultFile.write("\n%s\t"  % ExpeName)
	SNResultFile.write('\t'.join([str(MyScenario.Parameter(P, Silent=True))
								  for P in MyScenario.RelevantParamNames()]))
	SNResultFile.write("\t%s"  % LastStep)
	#SNResultFile.write("\t%d" % Step + 1)
	SNResultFile.write('\n' + "\t".join(Pop.Dump('SignalInvestment')))	  
	SNResultFile.write('\n' + "\t".join(Pop.Dump('DistanceToBestFriend')))	  
	SNResultFile.close()


		
def Start(BatchMode=False):
	if MyScenario.Parameter('BatchMode'):
		ExperienceID = strftime("%y%m%d%H%M%S")
	else:
		ExperienceID = strftime("_Interactive")
	Observer = SNetwork_Observer(ExperienceID)   # Observer contains statistics
	Pop = Population(Observer.NbAgents, Observer)   # set of flying swallows
	if BatchMode :
		for Step in range(MyScenario.Parameter('TimeLimit')):
			#print '.',
			Pop.One_Run()
			if os.path.exists('stop'):
				break
		Dump(Pop, Observer.ExperienceID, Observer.StepId)
		return

	# Interactive mode
	print 'RankEffects:', InterActions.RankEffects	
	" launching window "
	MainWindow = Tkinter.Tk()   # creates the main window
	MainWindow.title('Social Network')
	MainWindow.GParent = MainWindow

	Frame = SNetwork_Frame(MainWindow, Pop.One_Run, Observer)   

	# displaying the window
	MainWindow.lift()
	MainWindow.focus_force()
	# Entering main loop
	MainWindow.mainloop()   # control is given to the window system
	Dump(Pop, Observer.ExperienceID, Observer.StepId )
	


if __name__ == "__main__":

	BatchMode = MyScenario.Parameter('BatchMode')

	if BatchMode:
		Start(BatchMode=True)
	else:
		print __doc__
		Start(BatchMode=False)
		print "Bye......."
	
##      raw_input("\n[Return]")


__author__ = 'Dessalles'
