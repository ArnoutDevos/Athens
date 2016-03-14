#!/usr/bin/env python
################################################################################################
# Remplacement d'une chaine de carateres dans le contenu d'un fichier		 JLD 2012		 #
################################################################################################
"""	 Remplacement d'une chaine de carateres dans un fichier 
	ou dans les fichiers d'une arborescence (mais voir substree pour ca)
"""

import sys
import os.path
import time

def usage(cmd):
		print('Usage: %s <NomFich> <StringToFind> <NewString>' % cmd)
		time.sleep(2)
		
EscapeChar= {'\\t':'\t', '\\r':'\r', '\\n':'\n'}

def escapes(StrIn):
	if StrIn in EscapeChar:
		return EscapeChar[StrIn]
	return StrIn

def SubstituteInFile(FileName, OldString, NewString, Verbose=1, CommentLineChars=''):
	"""	Replaces OldString by NewString anywhere in File named FileName 
		except in lines starting with characters present in CommentLineChars
	"""
	try:
		Content = open(FileName).read()
		ContentLines = open(FileName).readlines()
	except IOError:
		print('Unable to open file {0}'.format(FileName))
		#raise SystemExit
		raise Exception
	if Content.find(escapes(OldString)) >= 0:
		# open(FileName,'w').write(Content.replace(escapes(OldString),escapes(NewString)))
		NewLines = []
		for Line in ContentLines:
			if Line[0] not in CommentLineChars:
				NewLines.append(Line.replace(escapes(OldString),escapes(NewString)))
			else:
				NewLines.append(Line)
		open(FileName,'w').write(''.join(NewLines))
		
		if Verbose:	print('Replacement done in {0}'.format(FileName))
	elif Verbose > 1:	print('String [{0}] not found in {1}'.format(OldString, FileName))


if __name__ == '__main__':

	print(__doc__)
	if len(sys.argv) > 3 and os.path.exists(sys.argv[-3]):	# ???
		SubstituteInFile(sys.argv[-3], sys.argv[-2], sys.argv[-1])
	else:
		usage(os.path.basename(sys.argv[0]))

	# Detabify
	# SubstituteInTree('.', 'py', '\\t', '    ', Verbose=False, CommentLineChars='#')
	# Tabify
	# SubstituteInTree('.', 'py', '    ', '\\t', Verbose=False, CommentLineChars='#')
			

__author__ = 'Dessalles'
