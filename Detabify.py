#!/usr/bin/env python
################################################################################################
# Remplacement d'une chaine de carateres dans le contenu d'un fichier		 JLD 2012		 #
################################################################################################
"""	 Replaces tabs by chars in python files
"""

import sys
import os.path
import time

sys.path.append('Tools')

import Walk

if __name__ == '__main__':
	print __doc__
	print("Do you want to replace all tabs by spaces in python source files")
	if raw_input('? ').lower().startswith('y'):
		# Detabify
		Walk.SubstituteInTree('.', '.*.py$', '\\t', ' '*4, CommentLineChars='#')
		# ReTabify
		# Walk.SubstituteInTree('.', '.*.py$', ' '*4, '\\t', CommentLineChars='#')
		print('Done')
	else:
		print ('Nothing done')

__author__ = 'Dessalles'
