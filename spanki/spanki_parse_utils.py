#!/usr/bin/env python
# encoding: utf-8

import re
import sys
import argparse
import pysam
import collections
import math


class Junctionid:
	"""
	Base class for a junction, from juncid
	Expects junction id
	"""
	def __init__(self, juncid):
		chr = juncid.split(':')[0]
		coords = juncid.split(':')[1]
		strand = juncid.split(':')[2]
		#start = int(coords.split('_')[0])
		start = int(coords.split('_')[0]) - 1
		end = int(coords.split('_')[1])
		self.chr = chr
		self.start = start
		self.end = end
		self.strand = strand.strip()
		self.intronsize = end - start
		#self.accid = str(coords.split('_')[0])
		if self.strand == "+":
			self.donid = chr + ":" + str(start + 1)
			self.donor = start
			self.accid = chr + ":" + str(end)
			self.donor = start
			self.acceptor = end
		elif self.strand == "-":
			self.accid = chr + ":" + str(start + 1)
			self.donid = chr + ":" + str(end)
			self.donor = end
			self.acceptor = start
		else:
			print "strand character is", self.strand
			quit("Don't recognize strand")
	def display(self):
		print self.intronsize
		print self.chr
		print self.donid
		#print self.dastring
		#print "start", self.start
		#print "L read anchor   ", self.left_read_anchor
		#print "L genome aligned", self.left_genome_aligned
		#print "R read anchor   ", self.right_read_anchor
		#print "R genome aligned", self.right_genome_aligned
		#print "READ            ", self.readseq


def parseRefAsBam(samfile):
	EDG = {}
	JUN = {}
	for alignedread in samfile:
		if (len(alignedread.cigar) > 1):
			strand = alignedread.tags[0][1] 
			mytid = alignedread.tid
			txid = str(alignedread.qname)
			chr = samfile.getrname(mytid)
			start = alignedread.pos
			offset = start
			gaps = []
			matches = []
			addnextmatch = 0
			addtomatch = 0
			for i in alignedread.cigar:
				# Iterate over each entry in cigar
				# matches should be appended 
				if (i[0] == 0):
					if addnextmatch > 0:
						matches[-1] = matches[-1] + i[1] + addtomatch
						addnextmatch = 0
						addtomatch = 0
					else: matches.append(i[1])
				elif (i[0] == 3):
					gaps.append(i[1])
				elif (i[0] == 1):
					addnextmatch += 1
					addtomatch = 0
				elif (i[0] == 2):
					addnextmatch += 1
					addtomatch = i[1]
									
			if len(gaps) > 0:
				for i in range(len(gaps)):
					junc_left_side = start + matches[i] + 1
					junc_right_side = start + matches[i] + gaps[i] 
					juncid = chr + ":" + str(junc_left_side) + "_" + str(junc_right_side) + ":" + strand
					if (strand == '+'):
						donid = chr + ":" + str(junc_left_side)
						accid = chr + ":" + str(junc_right_side)
					else:
						donid = chr + ":" + str(junc_right_side)
						accid = chr + ":" + str(junc_left_side)
					EDG[donid] = txid
					EDG[accid] = txid
					JUN[juncid] = 1
					start = start + matches[i] + gaps[i];	
					
 	#samfile.close()
 	refjunclist = JUN.keys()
 	refjunclist.sort()
 	return EDG, refjunclist


def edgeConnections(samfile):
	EDGconnections = collections.defaultdict(list)
	donlist = []
	acclist = []
	for alignedread in samfile:
		if (len(alignedread.cigar) > 1):
			strand = alignedread.tags[0][1] 
			mytid = alignedread.tid
			txid = str(alignedread.qname)
			chr = samfile.getrname(mytid)
			start = alignedread.pos
			offset = start
			gaps = []
			matches = []
			addnextmatch = 0
			addtomatch = 0
			for i in alignedread.cigar:
				# Iterate over each entry in cigar
				# matches should be appended 
				if (i[0] == 0):
					if addnextmatch > 0:
						matches[-1] = matches[-1] + i[1] + addtomatch
						addnextmatch = 0
						addtomatch = 0
					else: matches.append(i[1])
				elif (i[0] == 3):
					gaps.append(i[1])
				elif (i[0] == 1):
					addnextmatch += 1
					addtomatch = 0
				elif (i[0] == 2):
					addnextmatch += 1
					addtomatch = i[1]
									
			if len(gaps) > 0:
				for i in range(len(gaps)):
					junc_left_side = start + matches[i] + 1
					junc_right_side = start + matches[i] + gaps[i] 
					juncid = chr + ":" + str(junc_left_side) + "_" + str(junc_right_side) + ":" + strand
					j1 = Junctionid(juncid)
					donid = j1.donid
					accid = j1.accid
					#if (strand == '+'):
					#	donid = chr + ":" + str(junc_left_side)
					#	accid = chr + ":" + str(junc_right_side)
					#else:
					#	donid = chr + ":" + str(junc_right_side)
					#	accid = chr + ":" + str(junc_left_side)
					EDGconnections[donid].append(accid)
					EDGconnections[accid].append(donid)
					donlist.append(donid)
					acclist.append(accid)
					start = start + matches[i] + gaps[i];	
					
 	#samfile.close()
 	donlist = list(set(donlist))
 	acclist = list(set(acclist))
 	return EDGconnections, donlist, acclist

def edgeConnectionsPrevious(samfile):
	EDGconnections = collections.defaultdict(list)
	donlist = []
	acclist = []
	for alignedread in samfile:
		if (len(alignedread.cigar) > 1):
			strand = alignedread.tags[0][1] 
			mytid = alignedread.tid
			txid = str(alignedread.qname)
			chr = samfile.getrname(mytid)
			start = alignedread.pos
			offset = start
			gaps = []
			matches = []
			addnextmatch = 0
			addtomatch = 0
			for i in alignedread.cigar:
				# Iterate over each entry in cigar
				# matches should be appended 
				if (i[0] == 0):
					if addnextmatch > 0:
						matches[-1] = matches[-1] + i[1] + addtomatch
						addnextmatch = 0
						addtomatch = 0
					else: matches.append(i[1])
				elif (i[0] == 3):
					gaps.append(i[1])
				elif (i[0] == 1):
					addnextmatch += 1
					addtomatch = 0
				elif (i[0] == 2):
					addnextmatch += 1
					addtomatch = i[1]
									
			if len(gaps) > 0:
				for i in range(len(gaps)):
					junc_left_side = start + matches[i] + 1
					junc_right_side = start + matches[i] + gaps[i] 
					juncid = chr + ":" + str(junc_left_side) + "_" + str(junc_right_side) + ":" + strand
					if (strand == '+'):
						donid = chr + ":" + str(junc_left_side)
						accid = chr + ":" + str(junc_right_side)
					else:
						donid = chr + ":" + str(junc_right_side)
						accid = chr + ":" + str(junc_left_side)
					EDGconnections[donid].append(accid)
					EDGconnections[accid].append(donid)
					donlist.append(donid)
					acclist.append(accid)
					start = start + matches[i] + gaps[i];	
					
 	#samfile.close()
 	donlist = list(set(donlist))
 	acclist = list(set(acclist))
 	return EDGconnections, donlist, acclist

def get_junc_positions(samfile):
	#FBtr0091780	16	chr3R	4029917	255	523M54N59M	*	0	0	*	*	XS:A:-
	TX = collections.defaultdict(lambda : collections.defaultdict(dict))
	#JUNCTIONS = collections.defaultdict(lambda : collections.defaultdict(dict))
	JUNCORDERS = collections.defaultdict(list)
	JUNCSTARTD = collections.defaultdict(list)
	JUNCENDD = collections.defaultdict(list)
	#EDG = {}
	#JUN = {}
	for alignedread in samfile:
		if (len(alignedread.cigar) > 1):
			strand = alignedread.tags[0][1] 
			mytid = alignedread.tid
			txid = str(alignedread.qname)
			chr = samfile.getrname(mytid)
			start = alignedread.pos
			end = alignedread.aend
			
			if strand == "+":
				txstart = start
				txend = end
			elif strand == "-":
				txstart = end
				txend = start
			else:
				quit("Don't recongnize strand")
			
					
			
			numjuncs = len(alignedread.cigar) - 1
			TX[txid]['chr'] = chr
			TX[txid]['start'] = start
			TX[txid]['end'] = end
			TX[txid]['strand'] = strand
			TX[txid]['numjuncs'] = numjuncs
			juncorder = 0
			
			offset = start
			gaps = []
			matches = []
			addnextmatch = 0
			addtomatch = 0
			for i in alignedread.cigar:
				# Iterate over each entry in cigar
				# matches should be appended 
				if (i[0] == 0):
					if addnextmatch > 0:
						matches[-1] = matches[-1] + i[1] + addtomatch
						addnextmatch = 0
						addtomatch = 0
					else: matches.append(i[1])
				elif (i[0] == 3):
					gaps.append(i[1])
				elif (i[0] == 1):
					addnextmatch += 1
					addtomatch = 0
				elif (i[0] == 2):
					addnextmatch += 1
					addtomatch = i[1]
									
			if len(gaps) > 0:
				for i in range(len(gaps)):
					junc_left_side = start + matches[i] + 1
					junc_right_side = start + matches[i] + gaps[i] 
					juncid = chr + ":" + str(junc_left_side) + "_" + str(junc_right_side) + ":" + strand
					if (strand == '+'):
						donid = chr + ":" + str(junc_left_side)
						accid = chr + ":" + str(junc_right_side)
					else:
						donid = chr + ":" + str(junc_right_side)
						accid = chr + ":" + str(junc_left_side)
					juncorder += 1
					
					#JUNCTIONS[juncid]['order'].append(juncorder)
					#print "Juncorder is ", juncorder
					if strand == "+":
						startd = junc_left_side - txstart
						endd = txend - junc_right_side
						##print junc_left_side
						##print txstart
						##print txend
						##print startd
						##print endd
					elif (strand == "-"):
						startd = txstart - junc_right_side
						endd = junc_right_side - txend
					else:
						quit("Don't recongnize strand")
										
					##print "bleh", endd
					JUNCORDERS[juncid].append(juncorder)
					JUNCSTARTD[juncid].append(startd)
					JUNCENDD[juncid].append(endd)
					start = start + matches[i] + gaps[i];	


					
 	#samfile.close()
 	#refjunclist = JUN.keys()
 	#refjunclist.sort()
 	return TX, JUNCORDERS, JUNCSTARTD, JUNCENDD


def parseRefAsBamReturnJuncs(samfile):
	EDG = {}
	#JUN = {}
	JUN = collections.defaultdict(lambda : collections.defaultdict(dict))
	for alignedread in samfile:
		if (len(alignedread.cigar) > 1):
			strand = alignedread.tags[0][1] 
			mytid = alignedread.tid
			txid = str(alignedread.qname)
			chr = samfile.getrname(mytid)
			start = alignedread.pos
			offset = start
			gaps = []
			matches = []
			addnextmatch = 0
			addtomatch = 0
			for i in alignedread.cigar:
				# Iterate over each entry in cigar
				# matches should be appended 
				if (i[0] == 0):
					if addnextmatch > 0:
						matches[-1] = matches[-1] + i[1] + addtomatch
						addnextmatch = 0
						addtomatch = 0
					else: matches.append(i[1])
				elif (i[0] == 3):
					gaps.append(i[1])
				elif (i[0] == 1):
					addnextmatch += 1
					addtomatch = 0
				elif (i[0] == 2):
					addnextmatch += 1
					addtomatch = i[1]
									
			if len(gaps) > 0:
				for i in range(len(gaps)):
					junc_left_side = start + matches[i] + 1
					junc_right_side = start + matches[i] + gaps[i] 
					juncid = chr + ":" + str(junc_left_side) + "_" + str(junc_right_side) + ":" + strand
					if (strand == '+'):
						donid = chr + ":" + str(junc_left_side)
						accid = chr + ":" + str(junc_right_side)
					else:
						donid = chr + ":" + str(junc_right_side)
						accid = chr + ":" + str(junc_left_side)
					EDG[donid] = txid
					EDG[accid] = txid
					if JUN[juncid]:
						JUN[juncid].append(txid)
					else:
						JUN[juncid] = [txid]
					start = start + matches[i] + gaps[i];	
					
 	#samfile.close()
 	refjunclist = JUN.keys()
 	refjunclist.sort()
 	return EDG, refjunclist, JUN

def parseRefSelectionAsBam(samfile,txlist):
	EDG = {}
	JUN = {}
	for alignedread in samfile:
		mytid = alignedread.tid
		txid = str(alignedread.qname)
		if txid in txlist:
			if (len(alignedread.cigar) > 1):
				strand = alignedread.tags[0][1] 
				mytid = alignedread.tid
				txid = str(alignedread.qname)
				chr = samfile.getrname(mytid)
				start = alignedread.pos
				offset = start
				gaps = []
				matches = []
				addnextmatch = 0
				addtomatch = 0
				for i in alignedread.cigar:
					# Iterate over each entry in cigar
					# matches should be appended 
					if (i[0] == 0):
						if addnextmatch > 0:
							matches[-1] = matches[-1] + i[1] + addtomatch
							addnextmatch = 0
							addtomatch = 0
						else: matches.append(i[1])
					elif (i[0] == 3):
						gaps.append(i[1])
					elif (i[0] == 1):
						addnextmatch += 1
						addtomatch = 0
					elif (i[0] == 2):
						addnextmatch += 1
						addtomatch = i[1]
										
				if len(gaps) > 0:
					for i in range(len(gaps)):
						junc_left_side = start + matches[i] + 1
						junc_right_side = start + matches[i] + gaps[i] 
						juncid = chr + ":" + str(junc_left_side) + "_" + str(junc_right_side) + ":" + strand
						if (strand == '+'):
							donid = chr + ":" + str(junc_left_side)
							accid = chr + ":" + str(junc_right_side)
						else:
							donid = chr + ":" + str(junc_right_side)
							accid = chr + ":" + str(junc_left_side)
						EDG[donid] = txid
						EDG[accid] = txid
						JUN[juncid] = 1
						start = start + matches[i] + gaps[i];	
					
 	#samfile.close()
 	refjunclist = JUN.keys()
 	refjunclist.sort()
 	return EDG, refjunclist


if __name__ == "__main__":
    sys.exit(main())
