#!/usr/bin/env python
##############################################################################
# EVOLIFE  www.dessalles.fr/Evolife                    Jean-Louis Dessalles  #
#            Telecom ParisTech  2014                       www.dessalles.fr  #
##############################################################################


##############################################################################
#  xy_game                                                   #
##############################################################################

	##############################################


import sys
if __name__ == '__main__':  sys.path.append('../..')  # for tests

######################################
# specific variables and functions   #
######################################

import random
from Evolife.Scenarii.Default_Scenario import Default_Scenario
from Evolife.Tools.Tools import percent
import random



class Scenario(Default_Scenario):
	
	


	def initialization(self):
		self.level_of_trust = 98.0
		self.percentage_x = 0  
		self.average = 50
		self.score = 0   
		self.choice = []
		self.scores = [0] *6
		self.maximum =0
		self.num_of_rounds = 0
		self.maximal_score =0

	def make_choice(self):
	# make_choice makes the choice for all the different players, and puts their choice (x & y) in a list.
		choice=[]
		x_in_list = 100 - self.level_of_trust
		y_in_list = self.level_of_trust
		my_list = ['x'] * int(x_in_list) + ['y'] * int(y_in_list)
		for each in range(6): #5 can be changed in the number of players as parameter
			choice.append(random.choice(my_list))
		self.choice = choice
		return choice

	def calculate_percentage_of_x(self):
		temp = self.make_choice()
		number_of_x = temp.count('x')
		num_of_players = len(self.make_choice())
		self.percentage_x = 100 / num_of_players * number_of_x
		return self.percentage_x	

	def update_level_of_trust(self):
		if self.choice.count('x') > 0:
			temp_level = self.level_of_trust - 0.1
		if self.choice.count('x') == 0:
			temp_level = self.level_of_trust		
		a = self.percentage_x
		if temp_level>100:
			self.level_of_trust = 95
			return self.level_of_trust
		else:
			self.level_of_trust = temp_level

			return self.level_of_trust

	def calculate_individual_score(self):
		temp = self.choice
		temp_list =[]
		Xs = temp.count('x')
		if Xs >0:
			for i,j in enumerate(temp):
				if j == 'x':
					temp_list.append(i)
			for each in temp_list:
				self.scores[each] += 1
		return self.scores

	def calculate_average_score(self):
		self.average = sum(self.scores)/6
		return self.average

	def calculate_maximum_score(self):
		self.maximum = max(self.scores)
		return self.maximum


	def display_(self):
		disp = [('black', 'self.level_of_trust')]		
		disp +=  [('white','self.maximum')] 
		disp +=  [('red','self.average')]
		#disp +=  [('blue','self.maximal_score')]
		return disp		

	def local_display(self,VariableID):
		if VariableID == 'self.level_of_trust':
			return self.level_of_trust
		elif VariableID == 'self.maximum':
			return self.percentage_x
		elif VariableID == 'self.average':
			return self.average
		#elif VariableID == 'self.maximal_score':
		#	return self.maximal_score

	def maximal_scoring(self):
		score = self.num_of_rounds *6
		self.maximal_score = score
		return self.maximal_score

	def start_game(self, members):

		self.calculate_percentage_of_x()
		self.calculate_individual_score()
		self.update_level_of_trust()
		self.num_of_rounds = self.num_of_rounds + 1
		self.calculate_average_score()
		self.calculate_maximum_score()
		self.maximal_scoring()

		if self.num_of_rounds == 1000:
			self.level_of_trust = 95 
			self.percentage_x = 0  
			self.average = 50
			self.score = 0   
			self.choice = []
			self.scores = [0] *6
			self.maximum =0
			self.num_of_rounds = 0


if __name__ == "__main__":
	Scenario()
	print(__doc__ + '\n')
	#raw_input('[Return]')
	
