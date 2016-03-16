#!/usr/bin/env python
##############################################################################
# EVOLIFE  www.dessalles.fr/Evolife                    Jean-Louis Dessalles  #
#            Telecom ParisTech  2014                       www.dessalles.fr  #
##############################################################################



##############################################################################
#  S_MyScenario : a scenario that does nothing ! (and to be customized)      #
##############################################################################

	#   If your scenario is 'XXX', copy this file to S_XXX.py.
	#   Indicate your name, date, context and abstract.
	#   Change the scenario name in the Evolife Configuration Editor
	#   for the [Run] button to execute S_XXX.py.
	#   You may have to edit the Evolife Configuration File (EvolifeConfigTree.xml)
	#   (exit from the Configuration Editor first!)
	#   to add new parameters. You may use you favorite editor (e.g. Emacs or Notepad++)
	#   or a specialized xml editor such as "Serna XML Editor" (available for free).
	#   Insert a section for your scenario by copy and modifying an existing scenario section.
	#   Then run the Evolife Configuration Editor again
		

""" EVOLIFE: Empty Scenario: This scenario does nothing ! It is supposed to be customized.

	Evolife scenarii may rewrite several functions listed here :
	(those marked with '+' are called from 'Group.py')
	(those marked with 'o' are called from 'Observer.py')
	
	COPY-PASTE the required functions from Default_Scenario.py

	- initialization(self): allows to define local variables
	- genemap(self):	initialises the genes on the gene map (see 'Genetic_map.py')
	- phenemap(self):   defines a list of phenotypic character names (see 'Phenotype.py')
	+ season(self, year):   makes periodic actions like resetting parameters
	+ behaviour(self, BestIndiv, AvgIndiv):   defines a behaviour to be displayed
	+ life_game(self, members): defines a round of interactions - calls the five following functions
		- start_game(self, members):	group-level initialization before starting interactions
			- prepare(self, indiv): individual initialization before starting interactions
		- interaction(self, Indiv, Partner):	defines a single interaction 
			- partner(self, Indiv, members):	select a partner among 'members' that will interact with 'Indiv'
		- end_game(self, members):  an occasion for a closing round after all interactions
		- evaluation(self, Indiv):  defines how the score of an individual is computed
		- lives(self, members): converts scores into life points
	+ couples(self, members): returns a list of couples for procreation (individuals may appear in several couples!)- Calls the following function:
		- parents(self, candidates):	selects two parents from a list of candidates (candidate = (indiv, NbOfPotentialChildren))
	+ new_agent(self, child, parents): initializes newborns
	+ remove_agent(self, agent): action to be performed when an agent dies
	+ update_positions(self, members, groupID):	assigns a position to agents
	o default_view(self): says which windows should be open at start up
	o legends(self): returns a string to be displayed at the bottom ot the Legend window.
	o display_(self):   says which statistics are displayed each year
	o local_display_(self):   allows to display locally defined values
	o wallpaper(self, Window):	if one wants to display different backgrounds in windows

						************

	You may pick those function from actual scenarios as well.

"""

	#=============================================================#
	#  HOW TO MODIFY A SCENARIO: read Default_Scenario.py		 #
	#=============================================================#

import sys
if __name__ == '__main__':  sys.path.append('../..')  # for tests


from Evolife.Scenarii.Default_Scenario import Default_Scenario

import numpy as np
import ahkab
from ahkab import ahkab, circuit, time_functions
import pylab

######################################
# specific variables and functions   #
######################################

class Scenario(Default_Scenario):

	######################################
	# All functions in Default_Scenario  #
	# can be overloaded				  #
	######################################


	def genemap(self):
		""" Defines the name of genes and their position on the DNA.
		Accepted syntax:
		['genename1', 'genename2',...]:   lengths and coding are retrieved from configuration
		[('genename1', 8), ('genename2', 4),...]:   numbers give lengths in bits; coding is retrieved from configuration
		[('genename1', 8, 'Weighted'), ('genename2', 4, 'Unweighted'),...]:	coding can be 'Weighted', 'Unweighted', 'Gray', 'NoCoding'
		"""
		#return [('W1',8, 'Weighted'), ('L1',8, 'Weighted'), ('R1',8, 'Weighted')] 
		return [('R1',8, 'Weighted'), ('R2',8, 'Weighted')]
		
		#Put here the declaration of W1 L1 R W2 L2
		# W = [100e-9 800e-9] [m]
		# L = [1e-6 3e-6] [m]
		# R = [100 1000] [ohm]

	def phenemap(self):
		""" Defines the set of non inheritable characteristics
		"""
		return ['Character']	# Elements in phenemap are integers between 0 and 100 and are initialized randomly

	def evaluation(self, Indiv):

		score = self.score_resistiveDivider(Indiv)
		
		Indiv.score(score, FlagSet=True)
		
	def score_resistiveDivider(self, Indiv):
		# Avoid shortcircuited resistances by making minimal resistance = 1
		R1=Indiv.gene_value('R1')+1
		R2=Indiv.gene_value('R2')+1
		
		# This circuit has a trivial solution: R1 = min(int)+1 and R2 = max(int)+1
		mycircuit = circuit.Circuit(title="Resistive Divider Circuit")
		
		# Reference ground
		gnd = mycircuit.get_ground_node()
		
		# This makes sure nodes are uniquely defined
		n1 = mycircuit.create_node('n1')
		n2 = mycircuit.create_node('n2')
		
		# Add the resistors to the schematic
		mycircuit.add_resistor("R1", 'n1', 'n2', value=R1)
		mycircuit.add_resistor("R2", 'n2', gnd, value=R2)

		# Add a voltage source 
		mycircuit.add_vsource("V1", n1="n1", n2=gnd, dc_value=1, ac_value=1)
		
		# Resistors are passives with a flat frequency characteristic so an operating point (op) analysis is enough
		op_analysis = ahkab.new_op()

		# Start the simulation and store the result in <result>
		result = ahkab.run(mycircuit, op_analysis)
		
		# The figure of merit to maximize is the voltage at the intermediate node between R1 and R2
		score = result['op']['VN2']*100

		return score
	
	def score_ota(self, Indiv):
		
		W=Indiv.gene_value('W1')+1
		L=Indiv.gene_value('L1')+1
		R=Indiv.gene_value('R1')+1
		
		#dumbscore = (W/L)*R

		## Couple with spice simulator and extract GBW
		# Create new circuit
		mycircuit = circuit.Circuit(title="Operational Transconductance Amplifier Circuit")
		
		# Reference ground
		gnd = mycircuit.get_ground_node()
		n1 = mycircuit.create_node('n1')
		n2 = mycircuit.create_node('n2')
		
		#print W
		#print R
		mycircuit.add_resistor("R1", 'n1', 'n2', value=W)
		mycircuit.add_resistor("R2", 'n2', gnd, value=R)

		mycircuit.add_vsource("V1", n1="n1", n2=gnd, dc_value=4, ac_value=1)
		
		#op_analysis = ahkab.new_op()
		ac_analysis = ahkab.new_ac(start=1e3, stop=1e5, points=3)
		
		#opa = ahkab.new_op()
		r = ahkab.run(mycircuit, ac_analysis)
		
		#print r['op'].results
		#print r['ac']['VN2']
		#print np.abs(r['ac']['VN2'][1])
		
		#dumbscore = np.abs(r['ac']['VN2'][1])*100
		
		#dumbscore = r['op']['VN2']*10
		
		#Indiv.score(dumbscore, FlagSet=True)
		
		# some stupid behaviour, to be replaced
		#if self.Parameter('Parameter1') == self.Parameter('Parameter2'):
			#Indiv.score(10)
		
	def display_(self):
		""" Defines what is to be displayed. It offers the possibility
			of plotting the evolution through time of the best score,
			the average score, and the average value of the
			various genes defined on the DNA.
			It should return a list of pairs (C,X)
			where C is the curve colour and X can be
			'best', 'average', or any gene name as defined by genemap
		"""
		return [('blue','average'),('white','best')]		
###############################
# Local Test                  #
###############################

if __name__ == "__main__":
	print(__doc__ + '\n')
	SB = Scenario()
	raw_input('[Return]')
	
