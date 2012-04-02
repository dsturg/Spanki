#!/usr/bin/env python
# encoding: utf-8

import re
import sys
import argparse
import pysam
import collections
import math


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
