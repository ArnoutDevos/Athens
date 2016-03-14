#!/usr/bin/env python
##############################################################################
# EVOLIFE  www.dessalles.fr/Evolife                    Jean-Louis Dessalles  #
#            Telecom ParisTech  2014                       www.dessalles.fr  #
##############################################################################


##############################################################################
#  Global definitions                                                        #
##############################################################################

""" EVOLIFE: Global constants and parameters
"""

import sys
import re

if __name__ == '__main__':  sys.path.append('../..')  # for tests

from Evolife.Tools.Tools import FileAnalysis, error


#########################################
# Loding global parameters             #
#########################################

# a function to check whether a string represents a positive or negative integer
isInZ = lambda x: x.isdigit() or (len(x)>1 and x[0]=='-' and x[1:].isdigit())

def Num(x):
	try:	return int(x)
	except ValueError:	return float(x)


class Parameters:
	""" class Parameters: includes all modifiable parameters
	"""

	def __init__(self, CfgFile='Evolife.evo'):
		self.Params = self.txt_to_cfg(CfgFile)	# dictionnary of (parameter name, value)
		self.relevant = set()  # list of parameters that are actually used
		# print self

	def txt_to_cfg(self,CfgTxtFile):
		""" retrieves a configuration from a text file
		"""
		try:
			# reads lines with following syntax:
			# [<Prefix/>*]<NameOfParameter> <ParameterValue> [<comments>]
			# Numerical parameters
			Numerical = FileAnalysis(CfgTxtFile, r"^([^#]\S*/)?(\w+)\s+(-?\d+(\.\d*)?)\s")	
			Numerical = [(V[1], Num(V[2])) for V in Numerical]
			# NonNumerical parameters
			# Alphabcal = FileAnalysis(CfgTxtFile, "^([^#]\S*/)?(\w+)\s+([^-0-9/]\S*).*$")	
			Alphabcal = FileAnalysis(CfgTxtFile, "^([^#]\S*/)?(\w+)\s+(\S*).*$")	
			Alphabcal = [(V[1],V[2]) for V in Alphabcal if not set(V[2]) <= set('-0123456789.') ]
			#cfg = dict([(V[1],int(V[2])) for V in R])
			cfg = dict(Numerical + Alphabcal)
##			if len(cfg) < len(Numerical + Alphabcal):
##				error("Evolife_Parameters: duplicated parameter", str([V for V in Numerical + Alphabcal if V not in cfg]))
			return cfg
		except IOError:
			error("Evolife_Parameters: Problem accessing configuration file", CfgTxtFile)
		return None

	def cfg_to_txt(self, CfgTxtFile):
		""" stores parameters into a text file
		"""
		Filout = open(CfgTxtFile, "w")
		Filout.write('\n'.join([p + '\t' + str(self.Params[p])
								for p in sorted(self.relevant)]))
		Filout.close()
		
	def Parameter(self,ParamName, Silent=False, Optional=False):
		try:
			p = self.Params[ParamName]
		except KeyError:
			if Optional:
				return None
			error("Evolife_Parameters: Attempt to reach undefined parameter: ", ParamName)
		if not Silent:
			self.relevant.add(ParamName)
		return p

	def Param(self, ParamName):	return self.Parameter(ParamName)
	
	def addParameter(self, Param, Value):
		" Adds a new parameter "
		self.Params[Param] = Value
		
	def ParamNames(self):
		return [P for P in self.Params if isInZ(str(self.Params[P]))]

	def RelevantParamNames(self):
		return sorted(self.relevant)

	def ParamValues(self):
		return self.Params.values()

	def Relevant(self, ParamName):
		return ParamName in self.relevant

	def None2Default(self, ParamName, Default=0):
		if self.Parameter(ParamName, Optional=True) is None:	return	Default
		return self.Parameter(ParamName)
		
	def __repr__(self):
		return '\n'.join(sorted([k+' =\t'+str(self.Params[k]) for k in self.Params]))




#################################
# Loading global parameters	 #
#################################

##Evolife_Parameters = Parameters('Evolife.evo')


	
	
if __name__ == "__main__":
	print(__doc__ + '\n')
	# print Defs.__doc__ + '\n'
	Evolife_Parameters = Parameters('../Evolife_.evo')
	print(Evolife_Parameters)
	raw_input('\n[Return]')


__author__ = 'Dessalles'
