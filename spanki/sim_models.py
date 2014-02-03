#!/usr/bin/env python
# encoding: utf-8
"""
sim_models

Models for generating simulated reads
"""


from __future__ import division 

import re
import sys
import argparse
import pysam
import collections
import math
import os
import csv
import pkgutil

from pyfasta import Fasta

# Biopython, for revcomp
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

from collections import defaultdict

from datetime import datetime, date

def getMMNUMmodel(model,mybp,customdir):

	mmprob = []

	if (model == "random"):
		mmprob = [30,20,15,10,5,5,5,5,2.5,2.5]
		#mmprob.extend([.25] * (mybp - 10))
		
	elif (model == "NIST"):
		'''
		Load model data
		'''
		myresource = str('data/' + model + '_mmtotals.txt')
		data = pkgutil.get_data(__name__, myresource)
		lines = data.split('\n')
		for line in lines:
			values = line.rstrip().split("\t")
			try:
				mmprob.append(int(values[1]))
			except:
				pass
		'''
		Alternatively, read line by line
		'''
		#d = os.path.dirname(sys.modules[__name__].__file__)
		#infile = open(os.path.join(d, 'data/mmtotals.txt'), 'r')
		#for line in infile:
		#	print line
		#	values = line.rstrip().split("\t")
		#	mmprob.append(int(values[1]))
		#infile.close()

	elif (model == "dm3"):
		myresource = str('data/' + model + '_mmtotals.txt')
		data = pkgutil.get_data(__name__, myresource)
		lines = data.split('\n')
		for line in lines:
			values = line.rstrip().split("\t")
			try:
				mmprob.append(int(values[1]))
			except:
				pass

	elif (model == "flyheads"):
		myresource = str('data/' + model + '_mmtotals.txt')
		data = pkgutil.get_data(__name__, myresource)
		lines = data.split('\n')
		for line in lines:
			values = line.rstrip().split("\t")
			try:
				mmprob.append(int(values[1]))
			except:
				pass

	elif (model == "custom"):
		myresource = str(customdir + '/mmtotals.txt')
		data = open(os.path.join('', myresource), 'rb').read()
		lines = data.split('\n')
		for line in lines:
			values = line.rstrip().split("\t")
			try:
				mmprob.append(int(values[1]))
			except:
				pass
	
	elif (model == "errorfree"):
		mmprob = 0
		
	return mmprob

def getMMPOSmodel(model,mybp,customdir):

	mmposprob = []
	
	if (model == "random"):
		#mmposprob = [25,20,15,10,10,5,5,5,2,2]
		#mmposprob.extend([.75] * (mybp - 10))

		# Simple weighted probability, where prob of a mismatch increases with increased length
		# First 10 bases are weighted the same:
		myseed = 10
		mmprob = [myseed] * 10
		for x in range(myseed,mybp):
			mmprob.append(myseed + x)
		mmposprob = mmprob

	elif (model == "NIST"):
		myresource = str('data/' + model + '_mmcounts.txt')
		data = pkgutil.get_data(__name__, myresource)
		lines = data.split('\n')
		for line in lines:
			values = line.split("\t")
			try:
				mmposprob.append(int(values[1]))
			except:
				pass
	elif (model == "dm3"):
		myresource = str('data/' + model + '_mmcounts.txt')
		data = pkgutil.get_data(__name__, myresource)
		lines = data.split('\n')
		for line in lines:
			values = line.rstrip().split("\t")
			try:
				mmposprob.append(int(values[1]))
			except:
				pass
	elif (model == "flyheads"):
		myresource = str('data/' + model + '_mmcounts.txt')
		data = pkgutil.get_data(__name__, myresource)
		lines = data.split('\n')
		for line in lines:
			values = line.rstrip().split("\t")
			try:
				mmposprob.append(int(values[1]))
			except:
				pass
	elif (model == "custom"):
		myresource = str(customdir + '/mmcounts.txt')
		data = open(os.path.join('', myresource), 'rb').read()
		lines = data.split('\n')
		for line in lines:
			values = line.rstrip().split("\t")
			try:
				mmposprob.append(int(values[1]))
			except:
				pass
	elif (model == "errorfree"):
		mmposprob = [0] * mybp

	return mmposprob

def getMMTYPEmodel(model,customdir):

	mmtypeprob = collections.defaultdict(lambda : collections.defaultdict(dict))

	if (model == "random"):
		# Equal probability of each
		mmtypeprob['A']['C'] = 25
		mmtypeprob['A']['G'] = 25
		mmtypeprob['A']['T'] = 25
		mmtypeprob['A']['N'] = 25
		mmtypeprob['C']['A'] = 25
		mmtypeprob['C']['G'] = 25
		mmtypeprob['C']['T'] = 25
		mmtypeprob['C']['N'] = 25
		mmtypeprob['G']['C'] = 25
		mmtypeprob['G']['A'] = 25
		mmtypeprob['G']['T'] = 25
		mmtypeprob['G']['N'] = 25
		mmtypeprob['T']['C'] = 25
		mmtypeprob['T']['G'] = 25
		mmtypeprob['T']['A'] = 25
		mmtypeprob['T']['N'] = 25
	elif (model == "errorfree"):
		# Equal probability of each
		mmtypeprob['A']['C'] = 25
		mmtypeprob['A']['G'] = 25
		mmtypeprob['A']['T'] = 25
		mmtypeprob['A']['N'] = 25
		mmtypeprob['C']['A'] = 25
		mmtypeprob['C']['G'] = 25
		mmtypeprob['C']['T'] = 25
		mmtypeprob['C']['N'] = 25
		mmtypeprob['G']['C'] = 25
		mmtypeprob['G']['A'] = 25
		mmtypeprob['G']['T'] = 25
		mmtypeprob['G']['N'] = 25
		mmtypeprob['T']['C'] = 25
		mmtypeprob['T']['G'] = 25
		mmtypeprob['T']['A'] = 25
		mmtypeprob['T']['N'] = 25
	elif (model == "NIST"):
		myresource = str('data/' + model + '_mmtypes.txt')
		data = pkgutil.get_data(__name__, myresource)
		lines = data.split('\n')
		for line in lines:
			values = line.rstrip().split("\t")
			bases = values[0].split(">")
			try:
				mmtypeprob[bases[0]][str(bases[1])] = values[1]
			except:
				pass
	elif (model == "dm3"):
		myresource = str('data/' + model + '_mmtypes.txt')
		data = pkgutil.get_data(__name__, myresource)
		lines = data.split('\n')
		for line in lines:
			values = line.rstrip().split("\t")
			bases = values[0].split(">")
			try:
				mmtypeprob[bases[0]][str(bases[1])] = values[1]
			except:
				pass
	elif (model == "flyheads"):
		myresource = str('data/' + model + '_mmtypes.txt')
		data = pkgutil.get_data(__name__, myresource)
		lines = data.split('\n')
		for line in lines:
			values = line.rstrip().split("\t")
			bases = values[0].split(">")
			try:
				mmtypeprob[bases[0]][str(bases[1])] = values[1]
			except:
				pass
	elif (model == "custom"):
		myresource = str(customdir + '/mmtypes.txt')
		data = open(os.path.join('', myresource), 'rb').read()
		lines = data.split('\n')
		for line in lines:
			values = line.rstrip().split("\t")
			bases = values[0].split(">")
			try:
				mmtypeprob[bases[0]][str(bases[1])] = values[1]
			except:
				pass

	return mmtypeprob

def getQUALmodel(model,mybp,customdir):

	quals = collections.defaultdict(lambda : collections.defaultdict(dict))
	qualpos = collections.defaultdict(lambda : collections.defaultdict(dict))

	qualstring = ""
	if (model == "random"):
		# Equal probability of each
		qualstring = "G" * mybp
	elif (model == "errorfree"):
		# Equal probability of each
		qualstring = "G" * mybp
	elif (model == "NIST"):
		myresource = str('data/' + model + '_qualitiescounts.txt')
		data = pkgutil.get_data(__name__, myresource)
		lines = data.split('\n')
		for line in lines:
			values = line.rstrip().split("\t")
			if len(values) > 1:
				quals[values[0]] = []
				quals[values[0]].extend(values[1:len(values)])
				totalquals = len(values) - 1
		mypos = 0
		while mypos < totalquals: # iterate over positions
			temp = 0
			for i in quals.keys(): # iterate over quality scores
				if int(quals[i][mypos]) > temp:
					qual = i
					temp = int(quals[i][mypos])
			#print mypos, qual 
			mypos += 1
			qualstring += str(qual)
	elif (model == "dm3"):
		myresource = str('data/' + model + '_qualitiescounts.txt')
		data = pkgutil.get_data(__name__, myresource)
		lines = data.split('\n')
		for line in lines:
			values = line.rstrip().split("\t")
			if len(values) > 1:
				quals[values[0]] = []
				quals[values[0]].extend(values[1:len(values)])
				totalquals = len(values) - 1
		mypos = 0
		while mypos < totalquals: # iterate over positions
			temp = 0
			for i in quals.keys(): # iterate over quality scores
				if int(quals[i][mypos]) > temp:
					qual = i
					temp = int(quals[i][mypos])
			#print mypos, qual 
			mypos += 1
			qualstring += str(qual)
	elif (model == "flyheads"):
		myresource = str('data/' + model + '_qualitiescounts.txt')
		data = pkgutil.get_data(__name__, myresource)
		lines = data.split('\n')
		for line in lines:
			values = line.rstrip().split("\t")
			if len(values) > 1:
				quals[values[0]] = []
				quals[values[0]].extend(values[1:len(values)])
				totalquals = len(values) - 1
		mypos = 0
		while mypos < totalquals: # iterate over positions
			temp = 0
			for i in quals.keys(): # iterate over quality scores
				if int(quals[i][mypos]) > temp:
					qual = i
					temp = int(quals[i][mypos])
			#print mypos, qual 
			mypos += 1
			qualstring += str(qual)
	elif (model == "custom"):
		myresource = str(customdir + '/qualitiescounts.txt')
		data = open(os.path.join('', myresource), 'rb').read()
		lines = data.split('\n')
		for line in lines:
			values = line.rstrip().split("\t")
			if len(values) > 1:
				quals[values[0]] = []
				quals[values[0]].extend(values[1:len(values)])
				totalquals = len(values) - 1
		mypos = 0
		while mypos < totalquals: # iterate over positions
			temp = 0
			for i in quals.keys(): # iterate over quality scores
				if int(quals[i][mypos]) > temp:
					qual = i
					temp = int(quals[i][mypos])
			#print mypos, qual 
			mypos += 1
			qualstring += str(qual)

	return qualstring[0:mybp]

def getMMQUALmodel(model,customdir):
	'''
	Return consensus quality score for mismatched bases
	'''
	qual = ""
	if (model == "random"):
		# Equal probability of each
		qual = "#"
	elif (model == "errorfree"):
		# Equal probability of each
		qual = "#"
	elif (model == "NIST"):
		myresource = str('data/' + model + '_qualitiescounts.txt')
		data = pkgutil.get_data(__name__, myresource)
		lines = data.split('\n')
		temp = 0
		for line in lines:
			values = line.rstrip().split("\t")
			bases = values[0].split(">")
			try:
				if int(values[1]) > temp:
					qual = values[0]
					temp = int(values[1])
			except:
				pass
	elif (model == "dm3"):
		myresource = str('data/' + model + '_qualitiescounts.txt')
		data = pkgutil.get_data(__name__, myresource)
		lines = data.split('\n')
		temp = 0
		for line in lines:
			values = line.rstrip().split("\t")
			bases = values[0].split(">")
			try:
				if int(values[1]) > temp:
					qual = values[0]
					temp = int(values[1])
			except:
				pass
	elif (model == "flyheads"):
		myresource = str('data/' + model + '_qualitiescounts.txt')
		data = pkgutil.get_data(__name__, myresource)
		lines = data.split('\n')
		temp = 0
		for line in lines:
			values = line.rstrip().split("\t")
			bases = values[0].split(">")
			try:
				if int(values[1]) > temp:
					qual = values[0]
					temp = int(values[1])
			except:
				pass
	elif (model == "custom"):
		myresource = str(customdir + '/qualitiescounts.txt')
		data = open(os.path.join('', myresource), 'rb').read()
		lines = data.split('\n')
		temp = 0
		for line in lines:
			values = line.rstrip().split("\t")
			bases = values[0].split(">")
			try:
				if int(values[1]) > temp:
					qual = values[0]
					temp = int(values[1])
			except:
				pass
	return qual

if __name__ == "__main__":
    sys.exit(main())

