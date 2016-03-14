#!/usr/bin/env python

import os

def scrute(ff, NroL, NroW):
	FF = open(ff, 'r')
	Lines = FF.readlines()
	FF.close()
	Words = Lines[NroL].split()
	try:
		return Words[NroW]
	except IndexError:
		print "Erreur fichier %s" % ff
		return ''

for F in os.listdir('.'):
##    if scrute(F, 1, 20) != '40':
	if scrute(F, 1, 13) != '5':
		print F, scrute(F, 1, 13)
	


__author__ = 'Dessalles'
