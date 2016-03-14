#!/usr/bin/env python
################################################################################################
# Remplacement d'une chaine de carateres dans le contenu d'un fichier		 JLD 2012		 #
################################################################################################
"""	 Removes CR chars introduced by MsWindows
"""

import sys
import os.path

sys.path.append('Tools')

import Walk

if __name__ == '__main__':
	print __doc__
	print("Do you want to remove all CR in python source files")
	if raw_input('? ').lower().startswith('y'):
		Walk.SubstituteInTree('.', '.*.py$', '\\r', '', Verbose=False)
		print('Done')
	else:
		print ('Nothing done')

__author__ = 'Dessalles'
