#!/usr/bin/env python
##############################################################################
# EVOLIFE  www.dessalles.fr/Evolife                    Jean-Louis Dessalles  #
#            Telecom ParisTech  2014                       www.dessalles.fr  #
##############################################################################





""" S_Zip : a scenario that simplifies or complexifies strings.
	In the simplification case, the idea is to start from a random string,
	and then to mutate it until one gets a simpler string, which means
	a string that can be better compressed.
"""
	#=============================================================#
	#  HOW TO MODIFY A SCENARIO: read Default_Scenario.py		 #
	#=============================================================#

import zlib
import bz2

import sys
if __name__ == '__main__':  sys.path.append('../..')  # for tests

from Evolife.Scenarii.Default_Scenario import Default_Scenario

######################################
# specific variables and functions   #
######################################

class Scenario(Default_Scenario):
	""" S_Zip : a scenario that simplifies or complexifies strings.
		In the simplification case, the idea is to start from a random string,
		and then to mutate it until one gets a simpler string, which means
		a string that can be better compressed.
	"""

	######################################
	# Most functions below overload some #
	# functions of Default_Scenario	  #
	######################################

	def genemap(self):
		""" Defines the name of genes and their position on the DNA.
		Accepted syntax:
		['genename1', 'genename2',...]:   lengths and coding are retrieved from configuration
		[('genename1', 8), ('genename2', 4),...]:   numbers give lengths in bits; coding is retrieved from configuration
		[('genename1', 8, 'Weighted'), ('genename2', 4, 'Unweighted'),...]:	coding can be 'Weighted', 'Unweighted', 'Gray', 'NoCoding'
		"""
		return ['String']   # length is retrieved from configuration file

	def evaluation(self, Indiv):
		" defines how the score of an individual is computed "
		if Indiv.score() > 0:
			# the individual has already been evaluated
			return
		if self.Parameter('BitString'):
			# DNA is translated into a binary string
			StrDNA = ''
			for ii in range(0,self.Parameter('GeneLength'),8):
				StrDNA += chr(Indiv.read_DNA(ii,ii+8, coding = self.Parameter('Weighted')))
		else:
			# DNA is translated into a string of charaters '0' and '1'
			StrDNA = ''.join([str(b) for b in Indiv.get_DNA()])
		if self.Parameter('bz2'):
			compressor = bz2.compress
		else:
			compressor = zlib.compress
		if self.Parameter('Simplify'):
			Indiv.score(self.Parameter('GeneLength') - len(compressor(StrDNA)), FlagSet=True)
		else:
			Indiv.score(len(compressor(StrDNA)), FlagSet=True)
		return

	def default_view(self):	return ['Genomes']
		
	def display_(self):
		""" Defines what is to be displayed. It offers the possibility
			of plotting the evolution through time of the best score,
			the average score, any locally defined value,
			and the average value of the various genes and phenes.
			It should return a list of pairs (C,X)
			where C is the curve colour and X can be
			'best', 'average', 'local', any gene name defined in genemap
			or any phene defined in phenemap
		"""
		return [('white','best'), ('blue2','average')]


###############################
# Local Test                  #
###############################

if __name__ == "__main__":
	print(__doc__ + '\n')
	SB = Scenario()
	raw_input('[Return]')
	
