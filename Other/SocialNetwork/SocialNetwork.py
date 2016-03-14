#!/usr/bin/env python
##############################################################################
# SocialNetwork   www.dessalles.fr/Evolife              Jean-Louis Dessalles #
#                 Telecom ParisTech  2014                www.dessalles.fr    #
##############################################################################


""" Study of the role of signalling in the emergence of social networks:
Individuals must find a compromise between efforts devoted to signalling
and the probability of attracting followers.
Links are supposed to be symmetrical.
"""


from time import sleep
from random import sample, choice

import sys
sys.path.append('../../..')

import Evolife.Scenarii.Parameters 
import SocialSimulation as SSim
import Evolife.Ecology.Learner as EL
from Evolife.Tools.Tools import boost, percent, LimitedMemory


# Global elements
class Global:
	def __init__(self):
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
		
	def Param(self, ParameterName):	return self.Parameters.Parameter(ParameterName)


Gbl = Global()

		
		
class Interactions:
	" A few functions used to rule interactions between agents "

	def __init__(self, FCompetence):
		self.FCompetence = FCompetence  # Competence function that computes the competence of individuals
		self.RankEffects = []   # table of decreasing investments in friendship
		self.RankEffect(0)
	
	def RankEffect(self, Rank):
		""" computes a decreasing coefficient depending on one's rank
			in another agent's address book.
		"""
		if self.RankEffects == []:
			# Initializing the table of social time given to friend
			# depending on friend's rank

			#	-----------------------------------------------------------
			#   T: total amount of available time
			#   tn: social time devoted to nth friend
			#   Tn: occupied social time with n friends
			#   T1 = t1
			#   Tn = Tn-1 + tn * (1 - Tn-1 / T)
			#   Tn = Tn-1 (1 - tn / T) + tn
			# T controls overlap:
			# T= 1 ==> all social time is crowded within constant time
			#	   much overlap, more friends does not decrease
			#	   each friend's share significantly
			# T= 100  ==> no overlap, social times add to each other,
			#	   shares diminish as the number of friends increases
			#	-----------------------------------------------------------
			
			RkEffect = Gbl.Param('RankEffect')/100.0
			if RkEffect == 0:
				RkEffect = 0.000001
			if Gbl.Param('SocialOverlap'):	T = 100 * RkEffect / Gbl.Param('SocialOverlap')
			else:	T = 10000.0
			for n in range(Gbl.Param('MaxFriends') + 2):
				tn = (RkEffect) ** (n+1)	# friend #0 gets time = RkEffect; friend #1 gets time = RkEffect**2;
				if n == 0:	Tn = RkEffect
				else:	Tn = self.RankEffects[n-1][0] * (1-tn/T)  + tn
				self.RankEffects.append((Tn,tn))
				
		if Rank >= 0:
			try:	return self.RankEffects[Rank]
			except IndexError:
				error('S_Signalling: RankEffect', str('Rank == %d' % Rank))
		else:	return (0,0)

	def NegotiateOffer(self,Agent,Partner):
		""" returns the ranks Agent and Partner are ready to assign to each other
			in their respective address book. Agent's rank recursively depends on
			Partner's attitude towards Agent.
		"""
		# non-recursive version
		MaxOffer = 100
		OldAgentOffer = OldPartnerOffer = (0,0) # (Rank, Offer)
		AgentOffer = PartnerOffer = (0, MaxOffer)	# AgentOffer = Agent's offer to Partner
		while (OldAgentOffer, OldPartnerOffer) != (AgentOffer,PartnerOffer):
			(OldAgentOffer, OldPartnerOffer) = (AgentOffer, PartnerOffer)
			PartnerOffer = self.Offer(Partner, AgentOffer, Apparent=True)
			if PartnerOffer[0] < 0: return (0,0)
			AgentOffer = self.Offer(Agent, PartnerOffer, Apparent=True)
			#print 'Negotiation2: %s offers %d and %s offers %d (at ranks %d, %d)' \
			#	 % (Agent.id,AgentOffer[1],Partner.id,PartnerOffer[1],AgentOffer[0], PartnerOffer[0])
			if AgentOffer[0] < 0:	return (0,0)
		return (AgentOffer[1], PartnerOffer[1])			

	def SocialOffer(self, Competence, PartnerRank, nbFriends):	   
		""" An agent's social offer depends on its alleged or real competence,
			on the rank it offers in its address book, and on the number of friends already
			present in the address book (as it may influence available time)
		"""
		if PartnerRank < 0:	return 0
		rankEffect = self.RankEffect(PartnerRank)[1]	# rank in address book matters
		sizeEffect = self.RankEffect(1 + nbFriends)[0] # size of address book matters
		return float(Competence * rankEffect) / sizeEffect

	def Offer(self, Agent, (PartnerRankOffer, PartnerSocialOffer), Apparent=True):
		""" Agent is going to make an offer to Partner, based on Partner's offer
		"""
		OfferedRank = Agent.accepts(PartnerSocialOffer)
		if Gbl.Param('SocialSymmetry') > 0 and OfferedRank >= 0:
			# Social symmetry supposes that friends put themselves at identical levels in their address book
			OfferedRank = max(PartnerRankOffer, OfferedRank) # worst of the two ranks
		SocialOffer = self.SocialOffer(self.FCompetence(Agent, Apparent), OfferedRank, Agent.nbFriends())
		#print Agent.id, Agent.Signal, OfferedRank, SocialOffer
		return (OfferedRank, SocialOffer)
		
	def groom(self, Indiv, Partner):
		""" The two individuals negotiate partnership.
			First they signal their competence.
			Then, they make a "social offer" based on the rank
			(in their "address book") they are ready to assign to each other.
			Lastly, each one independently decides to join the other or not.
			cf. Dunbar's "grooming hypothesis"
		"""

		# new interaction puts previous ones into question
		if Indiv.follows(Partner):	Indiv.end_friendship(Partner)	# symmetrical splitting up

		# Negotiation takes place
		(IndivOffer, PartnerOffer) =  self.NegotiateOffer(Indiv, Partner)

		# social links are established accordingly
		if IndivOffer == 0 or PartnerOffer == 0:
			# One doesn't care about the other
##            print "\nNo deal: %s(%d)->%d, %s(%d)->%d" % \
##                (Indiv.id, Indiv.Signal, IndivOffer, Partner.id, Partner.Signal, PartnerOffer)
			return # the deal is not made

		if not Indiv.get_friend(IndivOffer, Partner, PartnerOffer):
			print("***** Scenario Signalling: Negotiation not respected")
			print Indiv.id, 'was accepted by', Partner.id,
			print 'with offer-:', IndivOffer
			print Indiv.id, sorted(Indiv.friends.performances()),
			# print sorted(Indiv.followers.performances())
			print Partner.id, sorted(Partner.friends.performances()),
			# print sorted(Partner.followers.performances())
			error('S_Signalling', "Negotiation not respected")
			return # the deal is not made

		return
		
class Inter_Actions(Interactions):
	""" Inherits from the Interaction class with a particular competence """
	
	def __init__(self):
		Interactions.__init__(self, self.LCompetence)
	
	def LCompetence(self, Indiv, Apparent=False):
		if Gbl.Param('Transparency'):	return Indiv.Competence(Apparent=False)
		else:	return Indiv.Competence(Apparent)

#------------------------------------------		
#Interaction scenario
Gbl.InterActions = Inter_Actions()
#------------------------------------------		
	
	
class Individual(SSim.Social_Individual):
	"   class Individual: defines what an individual consists of "

	def __init__(self, IdNb, maxQuality=100):
		# Learnable features
		Features = {'SignalInvestment': 0 }   # propensity to signal one's quality
		self.Signal = 0
		self.BestSignal = 0	# best signal value in memory
		SSim.Social_Individual.__init__(self, IdNb, features=Features, maxQuality=maxQuality, parameters=Gbl)

	def Reset(self):		# called by Learner when born again
		SSim.Social_Individual.Reset(self)
		self.Points = 0
		self.Risk = 0
		self.update()

	def update(self, infancy=True):
		Colour = 'brown'
		if infancy and not self.adult():	Colour = 'pink'
		# self.Signal = percent(self.Features['SignalInvestment'] * self.Quality)
		#self.Signal = Gbl.InterActions.FCompetence(self, Apparent=True)
		self.Signal = self.Competence(Apparent=True)
		BR = self.bestRecord()
		if BR:	self.BestSignal = BR[0]['SignalInvestment'] * self.Quality / 100.0
		else:	self.BestSignal = 0
		# self.Position = (self.id, self.Features['SignalInvestment'], Colour)
		self.Position = (self.Quality, self.BestSignal+1, Colour, 8)	# 8 == size of blob in display
		if Gbl.Param('Links') and self.best_friend():
			self.Position += (self.best_friend().Quality, self.best_friend().BestSignal+1, 21, 1)
			
	def Competence(self, Apparent=False):
		" returns the actual quality of an individual or its displayed version "
		if Apparent:
			BC = Gbl.Param('BottomCompetence')
			Comp = percent(100-BC) * self.Quality + BC
			# Comp = BC + self.Quality
			VisibleCompetence = percent(Comp * self.Features['SignalInvestment'])	   
			return self.Limitate(VisibleCompetence, 0, 100)
		else:
			# return Comp
			return self.Quality

	def Interact(self, Partner):
		Gbl.InterActions.groom(self, Partner)

	def assessment(self):
		self.Points -= percent(Gbl.Param('SignallingCost') * self.Features['SignalInvestment'])
		# self.restore_symmetry() # Checks that self is its friends' friend
		for (Rank,Friend) in enumerate(self.Friends()):
			AgentSocialOffer = Gbl.InterActions.SocialOffer(Gbl.InterActions.FCompetence(self,
																				Apparent=False),
														Rank, self.nbFriends())
			Friend.Risk = percent(Friend.Risk * (100 - percent(Gbl.Param('CompetenceImpact')
																	   * AgentSocialOffer)))
				
	def wins(self, Points):
		"   stores a benefit	"
		self.Points += 100 - self.Risk
		# self.Points = Gbl.InterActions.Saturate(self.Points)
		SSim.Social_Individual.wins(self, self.Points)
	
	def __str__(self):
		return "%s[%0.1f]" % (self.id, self.Signal)
		
class Population(SSim.Social_Population):
	" defines the population of agents "

	def __init__(self, NbAgents, Observer):
		" creates a population of agents "
		SSim.Social_Population.__init__(self, Gbl, NbAgents, Observer, IndividualClass=Individual)
		

	def season_initialization(self):
		for agent in self.Pop:
			# agent.lessening_friendship()	# eroding past gurus performances
			if self.Param('EraseNetwork'):	agent.forgetAll()
			agent.Points = 0
			agent.Risk = 100	# maximum risk
				
	def Dump(self, Slot):
		""" Saving investment in signalling for each adult agent
			and then distance to best friend for each adult agent having a best friend
		"""
		if Slot == 'SignalInvestment':
			#D = [(agent.Quality, "%2.03f" % agent.Features['SignalInvestment']) for agent in self.Pop if agent.adult()]
			#D += [(agent.Quality, " ") for agent in self.Pop if not agent.adult()]
			D = [(agent.Quality, "%2.03f" % agent.Features['SignalInvestment']) for agent in self.Pop]
		if Slot == 'DistanceToBestFriend':
			D = [(agent.Quality, "%d" % abs(agent.best_friend().Quality - agent.Quality)) for agent in self.Pop if agent.adult() and agent.best_friend() is not None]
			D += [(agent.Quality, " ") for agent in self.Pop if agent.best_friend() == None or not agent.adult()]
		return [Slot] + [d[1] for d in sorted(D, key=lambda x: x[0])]
		

		
def Start(BatchMode=False):
	Observer = SSim.Social_Observer(Gbl.Parameters)   # Observer contains statistics
	Observer.setOutputDir('___Results')
	Observer.recordInfo('DefaultViews',	['Field', 'Network'])	# Evolife should start with that window open
	Pop = Population(Gbl.Param('NbAgents'), Observer)   # population of agents
	
	SSim.Start(Pop, Observer, BatchMode=BatchMode)
	


if __name__ == "__main__":
	boost()   # A technical trick that sometimes provides impressive speeding up with Python up to 2.6
	Start(BatchMode=Gbl.Param('BatchMode'))



__author__ = 'Dessalles'
