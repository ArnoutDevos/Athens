#!/usr/bin/env python
##############################################################################
# EVOLIFE  www.dessalles.fr/Evolife                    Jean-Louis Dessalles  #
#            Telecom ParisTech  2014                       www.dessalles.fr  #
##############################################################################



##############################################################################
#  S_CircuitDesign : a scenario that uses a Genetic Algorithm for Analog	 #
#				  	 Circuit optimization. 									 #
# 		 	 17/03/2016					                        Arnout Devos #
#            Telecom ParisTech  2016                     www.arnoutdevos.net #
##############################################################################	

""" EVOLIFE: CircuitDesign Scenario: a scenario that uses a Genetic Algorithm 
			 for Analog Circuit optimization.

packages to import:
ahkab : The Spice simulator used
numpy : To have the absolute value function
scipy : To use interpolation for the 3dB bandwidth calculation

WATCH OUT: the package for the spice simulator (ahkab) is not standard, one 
		   needs to install it from: https://ahkab.github.io/ahkab/

INFO: Resistive divider example provides an easy way to understand the code 
	  and is well documented. This resistive divider simulation can be executed 
	  by switching to the resistive divider functions in genemap(self) and 
	  evaluation(self, Indiv).
"""

import sys
if __name__ == '__main__':  sys.path.append('../..')  # for tests


from Evolife.Scenarii.Default_Scenario import Default_Scenario

import numpy as np
from scipy import interpolate
from ahkab import ahkab, circuit, time_functions

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
		return [('W1',8, 'Weighted'), ('L1',8, 'Weighted'), ('R1',8, 'Weighted')] # For OTA simulation
		#return [('R1',8, 'Weighted'), ('R2',8, 'Weighted')] #For resistive divider simulation

	def phenemap(self):
		""" Defines the set of non inheritable characteristics
		"""
		return ['Character']	# Elements in phenemap are integers between 0 and 100 and are initialized randomly

	def evaluation(self, Indiv):

		#score = self.score_resistiveDivider(Indiv) # Resistive divider simulation
		score = self.score_ota(Indiv) # OTA simulation
		
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
		
		# Use physical values for Width, Length and Resistance in Si Units
		minimalW = 1e-6
		minimalL = 10e-9
		minimalR = 10
		
		W=(Indiv.gene_value('W1')+1)*minimalW
		L=(Indiv.gene_value('L1')+1)*minimalL
		R=(Indiv.gene_value('R1')+1)*minimalR

		## Couple with spice simulator and extract GBW
		# Create new circuit
		ota = circuit.Circuit(title="Operational Transconductance Amplifier Circuit")
		
		# Reference ground
		gnd = ota.get_ground_node()
		
		# This makes sure nodes are uniquely defined
		n1 = ota.create_node('n1') # Top of current source
		
		n2 = ota.create_node('n2') # Output-, Drain of MOS1 (left)
		n3 = ota.create_node('n3') # Output+, Drain of MOS2 (right)
		
		n4 = ota.create_node('n4') # Power supply
		
		n5 = ota.create_node('n5') # Bias of gate +
		n6 = ota.create_node('n6') # Bias of gate -
		
		print 'variables'
		print W
		print L
		print R
		
		# Add all the components to the circuit
		ota.add_resistor("R1", 'n4', 'n2', value=R)
		ota.add_resistor("R2", 'n4', 'n3', value=R)
		
		ota.add_capacitor('c1', n1='n2', n2=gnd, value=5e-12)
		ota.add_capacitor('c2', n1='n3', n2=gnd, value=5e-12)
		
		ota.add_model('ekv', 'nmos', dict(TYPE='n', VTO=.4, KP=10e-6))
		
		ota.add_mos('m1', nd='n2', ng='n5', ns='n1', nb=gnd, model_label='nmos', w=W, l=L)
		ota.add_mos('m2', nd='n3', ng='n6', ns='n1', nb=gnd, model_label='nmos', w=W, l=L)
		
		# Use a small signal amplitude
		Vac = 0.02
		ota.add_vsource("V1", n1="n5", n2=gnd, dc_value=0.5, ac_value=Vac)
		ota.add_vsource("V2", n1="n6", n2=gnd, dc_value=0.5, ac_value=-Vac)
		
		ota.add_vsource("V3", n1="n4", n2=gnd, dc_value=2, ac_value=0)
		
		ota.add_isource('ib', n1='n1', n2=gnd, dc_value=1e-3)
		
		# The circuit is frequency dependent, so an AC analysis is required. Bandwith will be a FOM.
		ac_analysis = ahkab.new_ac(start=1e6, stop=1e10, points=10)

		r = ahkab.run(ota, ac_analysis)
		
		# The output voltage between the two drains of the transistors is used.
		Output = np.abs(r['ac']['VN3']-r['ac']['VN2'])
		
		## Gain calculation
		Gain = Output[2]/Vac
		score = Gain
		
		
		## Bandwidth calculation
		# Normalize the output to the low frequency value and convert to array
		norm_out = np.abs(Output)/np.abs(Output).max()
		
		frequencies = r['ac']['f']
		
		print norm_out
		print frequencies
		
		# Find the 3dB frequency. This is defined as the Bandwidth
		x1 = norm_out
		y1 = frequencies
		
		# Combine lists into list of tuples
		points = zip(x1, y1)

		# Sort list of tuples by x-value
		points = sorted(points, key=lambda point: point[0])

		# Split list of tuples into two list of x values any y values
		x1, y1 = zip(*points)
		
		frequencies_interpolated = interpolate.interp1d(x1, y1)
		
		if 0.5 > min(x1):
			print frequencies_interpolated(0.5)
			Bandwidth = frequencies_interpolated(0.5)
		else:
			Bandwith = 0
		
		# The total score function is based on both Gain and Bandwidth
		score = Gain*Bandwidth
		print score
		
		return score
		
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
	
