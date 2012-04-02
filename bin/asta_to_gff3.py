#!/usr/bin/env python
# encoding: utf-8
"""
spankysplice.py

Input junction coverage, astalavista splicing definitions
Output compiled data by event
"""
from __future__ import division 

import operator # For sorting dictionaries

import re
import sys
import argparse
import pysam
import collections
import math
import numpy
import csv
import os

from pyfasta import Fasta

# Custom modules to import:
import spanky.spanky_parse_utils as spanky_parse_utils
import spanky.spanky_utils as spanky_utils

from datetime import datetime, date




def GTFtoDict(infn):
    #logging.info("Parsing the GTF file %s." %(infn))
    gtflines = []
    infile = open(infn, 'r')
    numlines = 0
    for line in infile:
        numlines = numlines + 1
        gtflines.append(parseGTFlineToDict(line))

	#print "Processed %d lines in %s."
    #logging.info("Processed %d lines in %s." %(numlines, infn))
    infile.close()

    return gtflines

def addAttributesToGTFline(linedict):
    attributes = linedict["attribute"].split(";")
    del attributes[-1]
    attributes = [x.strip() for x in attributes]
    
    for attrib in attributes:
        attr = attrib.strip().split(" ")
        linedict[attr[0]] = attr[1].strip("\"")
    return linedict

def parseGTFlineToDict(line):
    values = line.split("\t")
    keys = ["seqname", "source", "feature", "start", "end", "score",
            "strand", "unknown", "attribute"]
    linedict = dict(zip(keys, values))
    linedict = addAttributesToGTFline(linedict)
    return linedict



def crossTabulate(x):
	D = collections.defaultdict(int)
	for y in x:
		D[y] += 1
	x = D.keys()
	x.sort()
	for x in D.keys():
		print "\t", x, "\t", D[x]

def readCtab(fname):
	lines = csv.DictReader(open(fname, 'rb'), delimiter='\t')
	CTAB = {}
	for line in lines:
		CTAB[line['tracking_id']] = line
	return CTAB

def readJtab(jfile):
	"""
	Reads in junction table
	Assumes first line has fields
	"""
	JTAB = {}
	lines = csv.reader(open(jfile, 'rb'), delimiter='\t')
	linecount = 0
	for line in lines:
		if (linecount < 1):
			#keys = line.split(" ")
			keys = line
			#print keys
		else: 
			values = line
			linedict = dict(zip(keys, values))
			juncid = str(linedict['juncid']).rstrip(" ")
			JTAB[juncid] = linedict 
		linecount += 1
	return JTAB


def chain_to_exons(flankstring,chainstring,chr,strand,eventcode):
	
	chains = chainstring.split(",");  # List of sites
	flanks = flankstring.split(",");  # The flanks (Always 2?)
	# Skipped exon example:
	#flanks "110877^,111907-"; structure "0,1-2^"; splice_chain ",111005-111117^"

	if eventcode == "AltFE":
		temp = chains[0]
		pattern = re.compile('\d+\^')
		m = pattern.search(chains[0])
		try:
			donor1 = m.group(0)
		except:
			donor1 = "none"
		m = pattern.search(chains[1])
		try:
			donor2 = m.group(0)
		except:
			donor2 = "none"
		path1 = str(donor1) + str(flanks[1])
		path2 = str(donor2) + str(flanks[1])
	elif eventcode == "AltLE":
		pattern = re.compile('\d+\-')
		m = pattern.search(chains[0])
		try:
			acceptor1 = m.group(0)
		except:
			acceptor1 = "none"
		m = pattern.search(chains[1])
		try:
			acceptor2 = m.group(0)
		except:
			acceptor2 = "none"
		path1 = flanks[0] + acceptor1
		path2 = flanks[0] + acceptor2
	else:
		path1 = flanks[0] + chains[0] + flanks[1];
		path2 = flanks[0] + chains[1] + flanks[1];
	
	
	print "PATH1", eventcode, path1
	print "PATH2", eventcode, path2
	# Traverse through each path.  Make a join for all cases where there is [0-9]-[0-9]^
	# Add correction to standardized id

	p = re.compile('\d+\^\d+\-')
	joins1 = []
	coords = p.findall(path1)
	for coord in coords:
		q = re.compile('\d+')
		x = q.findall(coord)
		if x[0] < x[1]: joins1.append(chr + ":" + str(int(x[0]) + 1) + "_" + str(int(x[1]) - 1) + ":" + strand)
		else: joins1.append(chr + ":" + str(int(x[1]) + 1) + "_" + str(int(x[0]) - 1) + ":" + strand)
	joins2 = []
	coords = p.findall(path2)
	for coord in coords:
		q = re.compile('\d+')
		x = q.findall(coord)
		if x[0] < x[1]: joins2.append(chr + ":" + str(int(x[0]) + 1) + "_" + str(int(x[1]) - 1) + ":" + strand)
		else: joins2.append(chr + ":" + str(int(x[1]) + 1) + "_" + str(int(x[0]) - 1) + ":" + strand)
		
	if len(joins1) < 1: 
		if flanks[0] == "null":
			joins1.append("none")
		else:
			joins1.append("irt")
	if len(joins2) < 1: joins2.append("none")
	joinstring = ",".join(joins1) + ";" + ",".join(joins2) 

	#Retained intron
	#chr2L	Undefined	as_event	155333	155784	1.39999	+	.	transcript_id "FBtr0078118,FBtr0301887"; gene_id "chr2L:155333-157666W"; flanks "155333(,155784^"; structure "0,1^2-"; splice_chain ",155429^155546-"; sources "Undefined,Undefined"; dimension "2_5"; degree "4"; localization "5UTR-CDS_maxCDS"; 
	#chr2L	Undefined	as_event	155333	155784	1.39999	+	.	transcript_id "FBtr0078118,FBtr0301886"; gene_id "chr2L:155333-157666W"; flanks "155333(,155784^"; structure "0,1^2-"; splice_chain ",155410^155466-"; sources "Undefined,Undefined"; dimension "2_5"; degree "4"; localization "5UTR-CDS_maxCDS"; 
	#chrX	Undefined	as_event	20344594	20345509	2.00000	-	.	transcript_id "FBtr0077320,FBtr0301537"; gene_id "chrX:20342674-20362133C"; flanks "20345509-,20344594^"; structure "0,1^2-"; splice_chain ",20345249^20345035-"; sources "Undefined,Undefined"; dimension "2_2"; degree "4"; localization "CDS_maxCDS"; 
	
	#Note:  None of them have "null" in the flanks field
	
	#Exon skip
	#chr2L	Undefined	as_event	110877	111907	1.13333	+	.	transcript_id "FBtr00                 "; gene_id "chr2L:106903-114433W"; flanks "110877^,111907-"; structure "0,1-2^"; splice_chain ",111005-111117^"; sources "Undefined,Undefined"; dimension "2_2"; degree "4"; localization "CDS_maxCDS"; 
	
	#No junction to assay
	#chrX	Undefined	as_event	21484626	21486039	2.00000	-	.	transcript_id "FBtr0070053,FBtr0070054"; gene_id "chrX:21482818-21486039C"; flanks "null,21484626^"; structure "1(2^4-,3("; splice_chain "21486039(21485651^21484978-,21485062("; sources "Undefined,Undefined"; dimension "2_2"; degree "4"; localization "5UTR-CDS_maxCDS"; 
	
	return joinstring;

def chain_to_joins(flankstring,chainstring,chr,strand,eventcode):
	
	chains = chainstring.split(",");  # List of sites
	flanks = flankstring.split(",");  # The flanks (Always 2?)
	# Skipped exon example:
	#flanks "110877^,111907-"; structure "0,1-2^"; splice_chain ",111005-111117^"

	if eventcode == "AltFE":
		temp = chains[0]
		pattern = re.compile('\d+\^')
		m = pattern.search(chains[0])
		try:
			donor1 = m.group(0)
		except:
			donor1 = "none"
		m = pattern.search(chains[1])
		try:
			donor2 = m.group(0)
		except:
			donor2 = "none"
		path1 = str(donor1) + str(flanks[1])
		path2 = str(donor2) + str(flanks[1])
	elif eventcode == "AltLE":
		pattern = re.compile('\d+\-')
		m = pattern.search(chains[0])
		try:
			acceptor1 = m.group(0)
		except:
			acceptor1 = "none"
		m = pattern.search(chains[1])
		try:
			acceptor2 = m.group(0)
		except:
			acceptor2 = "none"
		path1 = flanks[0] + acceptor1
		path2 = flanks[0] + acceptor2
	else:
		path1 = flanks[0] + chains[0] + flanks[1];
		path2 = flanks[0] + chains[1] + flanks[1];
	
	# Traverse through each path.  Make a join for all cases where there is [0-9]-[0-9]^
	# Add correction to standardized id

	p = re.compile('\d+\^\d+\-')
	joins1 = []
	coords = p.findall(path1)
	for coord in coords:
		q = re.compile('\d+')
		x = q.findall(coord)
		if x[0] < x[1]: joins1.append(chr + ":" + str(int(x[0]) + 1) + "_" + str(int(x[1]) - 1) + ":" + strand)
		else: joins1.append(chr + ":" + str(int(x[1]) + 1) + "_" + str(int(x[0]) - 1) + ":" + strand)
	joins2 = []
	coords = p.findall(path2)
	for coord in coords:
		q = re.compile('\d+')
		x = q.findall(coord)
		if x[0] < x[1]: joins2.append(chr + ":" + str(int(x[0]) + 1) + "_" + str(int(x[1]) - 1) + ":" + strand)
		else: joins2.append(chr + ":" + str(int(x[1]) + 1) + "_" + str(int(x[0]) - 1) + ":" + strand)
		
	if len(joins1) < 1: 
		if flanks[0] == "null":
			joins1.append("none")
		else:
			joins1.append("irt")
	if len(joins2) < 1: joins2.append("none")
	joinstring = ",".join(joins1) + ";" + ",".join(joins2) 

	#Retained intron
	#chr2L	Undefined	as_event	155333	155784	1.39999	+	.	transcript_id "FBtr0078118,FBtr0301887"; gene_id "chr2L:155333-157666W"; flanks "155333(,155784^"; structure "0,1^2-"; splice_chain ",155429^155546-"; sources "Undefined,Undefined"; dimension "2_5"; degree "4"; localization "5UTR-CDS_maxCDS"; 
	#chr2L	Undefined	as_event	155333	155784	1.39999	+	.	transcript_id "FBtr0078118,FBtr0301886"; gene_id "chr2L:155333-157666W"; flanks "155333(,155784^"; structure "0,1^2-"; splice_chain ",155410^155466-"; sources "Undefined,Undefined"; dimension "2_5"; degree "4"; localization "5UTR-CDS_maxCDS"; 
	#chrX	Undefined	as_event	20344594	20345509	2.00000	-	.	transcript_id "FBtr0077320,FBtr0301537"; gene_id "chrX:20342674-20362133C"; flanks "20345509-,20344594^"; structure "0,1^2-"; splice_chain ",20345249^20345035-"; sources "Undefined,Undefined"; dimension "2_2"; degree "4"; localization "CDS_maxCDS"; 
	
	#Note:  None of them have "null" in the flanks field
	
	#Exon skip
	#chr2L	Undefined	as_event	110877	111907	1.13333	+	.	transcript_id "FBtr00                 "; gene_id "chr2L:106903-114433W"; flanks "110877^,111907-"; structure "0,1-2^"; splice_chain ",111005-111117^"; sources "Undefined,Undefined"; dimension "2_2"; degree "4"; localization "CDS_maxCDS"; 
	
	#No junction to assay
	#chrX	Undefined	as_event	21484626	21486039	2.00000	-	.	transcript_id "FBtr0070053,FBtr0070054"; gene_id "chrX:21482818-21486039C"; flanks "null,21484626^"; structure "1(2^4-,3("; splice_chain "21486039(21485651^21484978-,21485062("; sources "Undefined,Undefined"; dimension "2_2"; degree "4"; localization "5UTR-CDS_maxCDS"; 
	
	return joinstring;


def geneidsFromTxlist(txlist,lookup):
	geneids = {}
	genenames = {}
	for tx in txlist:
		#print "** Checking tx ", tx
		geneids[lookup[tx]["gene_id"]] = 1
		genenames[lookup[tx]["gene_name"]] = 1
	gidlist = geneids.keys()
	gids = str(','.join(gidlist))
	gnamelist = genenames.keys()
	gnames = str(','.join(gnamelist))
	return gids,gnames


def summarize_events(eventdict,outfile):

	namedevents = collections.defaultdict(lambda : collections.defaultdict(dict))	
	unclassifiedevents = collections.defaultdict(lambda : collections.defaultdict(dict))	
	
	for event in eventdict.keys():
		"""
		Iterate over each item in eventdict
		"""	
		astacode = eventdict[event]['eventcode']
		if astacode in namedevents.keys():
			namedevents[astacode] += 1
		else:
			namedevents[astacode] = 1
		if astacode == "Unclassified":
			if eventdict[event]['structure'] in unclassifiedevents.keys():
				unclassifiedevents[eventdict[event]['structure']] += 1
			else:
				unclassifiedevents[eventdict[event]['structure']] = 1

	print  >> outfile, "Here is the breakdown for named events:"
	for code in namedevents.keys():
		print >> outfile, "\t", code, "-->", namedevents[code]
	print  >> outfile, "These are unclassified structures that have >= 10 events"
	sorted_x = sorted(unclassifiedevents.iteritems(), key=operator.itemgetter(1))
	totalnoname = 0
	for y in sorted_x:
		if (y[1] >= 10):
			print  >> outfile, "\t**", y[0], '\t', y[1]
			#quit()
		totalnoname += y[1]
	print >> outfile, "Total EVENTS with no english translation:  ", totalnoname
	print >> outfile, "Total STRUCTURES with no english translation:  ", len(unclassifiedevents.keys())

def countEventTypes(ECODES, astalines):

	Ks = ECODES.keys()
	
	EVENTCOUNT = {}
	NONAMEEVENTS = {}
	
	for line in astalines:
		"""
		Iterate over each line in asta def file
		"""	
		#print line['gene_id']
		astacode = line['structure']
		if astacode in Ks:
			eventcode = ECODES[astacode]
			#print astacode, "is", eventcode
			if EVENTCOUNT.has_key(eventcode): 
				EVENTCOUNT[eventcode] += 1
			else:
				EVENTCOUNT[eventcode] = 1
		else:
			if NONAMEEVENTS.has_key(astacode): 
				NONAMEEVENTS[astacode] += 1
			else:
				NONAMEEVENTS[astacode] = 1
	return EVENTCOUNT, NONAMEEVENTS

def print_events(NAMED,UNNAMED):
	print >> sys.stderr, "[%s] Counts of event types:" % (spanky_utils.timestamp())
	print >> sys.stderr, "[%s] Frequent occuring non-categorized:" % (spanky_utils.timestamp())
	sorted_x = sorted(UNNAMED.iteritems(), key=operator.itemgetter(1))
	totalnoname = 0
	totalevents = 0
	for y in sorted_x:
		if (y[1] >= 100):
			print y[0], '\t', y[1]
		totalnoname += y[1]
	print "Total non-characterized:  ", totalnoname
	for key in NAMED.keys():
		print key, '=>', NAMED[key]
		totalevents += NAMED[key]
 	totalevents += totalnoname
 	print "There are ", totalevents, " total events"

	# A sample entry
	#{'attribute': 'transcript_id "FBtr0089158,FBtr0089162"; gene_id "chr4:89956-129430W"; flanks "null,92947-"; structure "1(2^,3(4^"; splice_chain "89956(90020^,91032(91193^";
	#sources "Undefined,Undefined"; dimension "2_6"; degree "4"; localization "5UTR_maxCDS"; \n', 'joins2': ['chr4:91194_92946:+'], 'joins1': ['chr4:90021_92946:+'], 
	#'sources': 'Undefined,Undefined', 'splice_chain': '89956(90020^,91032(91193^', 'seqname': 'chr4', 'end': '92947', 'eventcode': 'altFE', 'source': 'Undefined', 'unknown': '.',
	#'feature': 'as_event', 'start': '89956', 'score': '1.22222', 'strand': '+', 'gnames': 'FBgn0085432', 'degree': '4', 'localization': '5UTR_maxCDS', 'gene_id': 'chr4:89956-129430W', 
	#'transcript_id': 'FBtr0089158,FBtr0089162', 'flanks': 'null,92947-', 'structure': '1(2^,3(4^', 'gids': 'pan', 'txlist2': ['FBtr0089162'], 'txlist1': ['FBtr0089158', 'FBtr0089162'], 
	#'joinstring': 'chr4:90021_92946:+;chr4:91194_92946:+', 'dimension': '2_6'}

		
def astacode_to_english(astacode):
	eventcode = ""
	fepattern = re.compile('\(')
	lepattern = re.compile('\)')
	ECODES = {}
	ECODES["1-2),3-4)"] = "AltLE"
	ECODES["1(2^,3(4^"] = "AltFE"
	ECODES["1-2^,3-4^"] = "mutexcl"
	ECODES["0,1-2^"] = "exonskip"
	ECODES["1^,2^"] = "altdonor"
	ECODES["1-,2-"] = "altacceptor"
	ECODES["0,1^2-"] = "retintron"
	ECODES["0,1-2^3-4^"] = "skip2exons"
	ECODES["1^4-,2^3-"] = "altdon_altacc"
	ECODES["1^3-,2^4-"] = "altdon_altacc"
	
	if astacode in ECODES.keys():
		eventcode = ECODES[astacode]
	else:
		fe = fepattern.search(astacode)
		le = lepattern.search(astacode)
		if (fe > 0): eventcode = "AltFE"
		elif (le > 0): eventcode = "AltLE"
		else: 
			eventcode = "Unclassified"
		if (fe > 0) & (le > 0): 
			eventcode = "FELE"
	
	return eventcode


def parse_asta_defs(line,lookup):
	"""
	Classify event types
	Identify inclusion/exclusion paths
	Generate "joinstrings"
	Add geneid using lookup table
	"""
	linedict = {}

	flankstring = line['flanks']
	chainstring = line['splice_chain']
	chr = line['seqname']
	astacode = line['structure']
	strand = line['strand']
	transcripts = line['transcript_id']
	astagene = line['gene_id']
	"""
	Get geneid from txids
	"""
	txie = transcripts.split(',')
	txlist1 = txie[0].split('/')
	txlist2 = txie[1].split('/')
	txlist = txlist1
	txlist.extend(txlist2)
	gnames, gids = geneidsFromTxlist(txlist,lookup)
	
	eventcode = astacode_to_english(astacode)
			
	joinstring = chain_to_joins(flankstring,chainstring,chr,strand,eventcode)

	"""
	Excluding some events
	(Be conservative)
	"""
	
	if (astacode == "0,0") or (eventcode == "FELE"):
		pass
	else:
		'''
		Instantiate events
		'''
		linedict = line
		linedict['joinstring'] = joinstring
		joins = joinstring.split(";")
		joinstring1 = joins[0]
		joinstring2 = joins[1]
		joins1 = joinstring1.split(',')
		joins2 = joinstring2.split(',')
		'''
		Adjustment for alt LE
		Only use most 5' join
		'''
		if eventcode == "AltLE":
			linedict['joins1'] = [joins1[0]]
			linedict['joins2'] = [joins2[0]]		
		else:
			linedict['joins1'] = joins1
			linedict['joins2'] = joins2
		linedict['txlist1'] = txlist1
		linedict['txlist2'] = txlist2
		linedict['eventcode'] = eventcode
		linedict['gnames'] = gnames
		linedict['gids'] = gids
	return linedict


def get_coverage(eventdict):
	#DSX joinid
	#testid = "chr3R:3761376_3761489:-"
	psi_errors = 0
	adj_ok = 0
	for eventid in eventdict.keys():	
		cov1 = 0
		cov2 = 0
		joins1 = eventdict[eventid]['joins1']
		joins2 = eventdict[eventid]['joins2']
		j1 = []
		j2 = []
		#print joins1, joins1

		join1COV = []
		join2COV = []
		sites1 = len(joins1)
		sites2 = len(joins2)

		for join in joins1:
			if join == 'irt':
				sites1 = sites2 * 2
				try:
					tempcov = 0
					for x in joins2:
						tempcov += int(JTAB[x]['irt'])
					join1COV.append(tempcov)
				except KeyError:
					join1COV.append(0)
			else:
				try:
					join1COV.append(int(JTAB[join]['cov']))
				except KeyError:
					join1COV.append(0)

		for join in joins2:
			try:
				join2COV.append(int(JTAB[join]['cov']))
			except KeyError:
				join2COV.append(0)
	
		'''
		# For applying a correction to multiexon interior events
		'''

		##################################
		# Alternate version
		# by subtracting neighbor coverage
		##################################
		#interior_events = ["Unclassified","mutexcl","exonskip","skip2exons"]
		#join1_adjust = 0
		#join2_adjust = 0
		#if eventdict[eventid]['eventcode'] in interior_events:
		#	if (len(joins1) > 1):
		#		for join in joins1:
		#			try:
		#				join1_adjust += int(JTAB[join]['ancov']) 
		#				join1_adjust += int(JTAB[join]['dncov'])
		#			except KeyError:
		#				pass
		#		try:
		#			join1_adjust -= int(JTAB[joins1[0]]['dncov']) 
		#		except:
		#			pass
		#		try:
		#			join1_adjust -= int(JTAB[joins1[-1]]['ancov'])
		#		except:
		#			pass
		#	if (len(joins2) > 1):
		#		for join in joins2:
		#			try:
		#				join2_adjust += int(JTAB[join]['ancov']) 
		#				join2_adjust += int(JTAB[join]['dncov'])
		#			except:
		#				pass
		#		try:
		#			if (eventdict[eventid]['strand'] == '+'): join2_adjust -= int(JTAB[joins2[0]]['dncov']) 
		#			else: join2_adjust -= int(JTAB[joins2[0]]['ancov']) 
		#		except:
		#			pass
		#		try:
		#			if (eventdict[eventid]['strand'] == '+'): join2_adjust -= int(JTAB[joins2[-1]]['ancov'])
		#			else: join2_adjust -= int(JTAB[joins2[-1]]['dncov']) 
		#		except:
		#			pass
					
		
		inc = sum(join1COV)
		exc = sum(join2COV)
		
		if (eventdict[eventid]['eventcode'] == "AltFE") or (eventdict[eventid]['eventcode'] == "AltFE"):
			inc_adj = join1COV[0]
			exc_adj = join2COV[0]
		else:
			inc_adj = inc
			exc_adj = exc
		
		# If you apply the correction by subtraction
		#inc_adj = sum(join1COV) - join1_adjust
		#exc_adj = sum(join2COV) - join2_adjust
		#if (sum(join1COV) >= join1_adjust) and (sum(join2COV) >= join2_adjust):
		#	adj_ok += 1
		#else:
		#	psi_errors += 1
		
		eventdict[eventid]['inc'] = inc
		eventdict[eventid]['exc'] = exc

		if inc_adj > 0: eventdict[eventid]['inc_adj'] = inc_adj
		else: eventdict[eventid]['inc_adj'] = 0

		if exc_adj > 0: eventdict[eventid]['exc_adj'] = exc_adj
		else: eventdict[eventid]['exc_adj'] = 0

		try:
			psi = (inc / sites1) / ((inc / sites1) + (exc / sites2))
		except ZeroDivisionError:
			psi = 0
		try:
			psi_adj = (inc_adj) / ((inc_adj) + (exc_adj))
		except ZeroDivisionError:
			psi_adj = 0


		eventdict[eventid]['psi'] = "%.3f" % psi
		eventdict[eventid]['psi_adj'] = "%.3f" % psi_adj
		
		tx1 = eventdict[eventid]['txlist1']
		tx2 = eventdict[eventid]['txlist2']
		fpkm1 = 0
		fpkm2 = 0
		for tx in tx1:
			try: 
				fpkm1 += float(CTAB[tx]['FPKM'])
			except KeyError:
				fpkm1 += 0
		for tx in tx2:
			try:
				fpkm2 += float(CTAB[tx]['FPKM'])
			except KeyError:
				fpkm2 += 0
		eventdict[eventid]['fpkm1'] = "%.3f" % fpkm1
		eventdict[eventid]['fpkm2'] = "%.3f" % fpkm2
		try:
			eventdict[eventid]['fpkm_proportion'] = "%.3f" % (fpkm1 / (fpkm1 + fpkm2))
		except:
			eventdict[eventid]['fpkm_proportion'] = "Low data"
	#print "\tErrors:", psi_errors
	#print "\tOK", adj_ok
	
def print_output_table(mydict,handle,fieldorder):
	print >> handle, "id\t", '\t'.join(map(str,fieldorder))
	for myid in sorted(mydict.keys()):
		vals = []
		for field in fieldorder:
			vals.append(mydict[myid][field])
		print >> handle, myid, "\t", '\t'.join(map(str,vals))


def uniq(seq): 
	# Not order preserving
	keys = {} 
	for e in seq:
		keys[e] = 1 
	return keys.keys()

def gtf_to_tx_edges_dict(infn):
	'''
	Returns a dict of the edges
	keyed by txid_edge1 -- edge2
	For building GFF3 files
	'''
	txedgedict = collections.defaultdict(lambda : collections.defaultdict(dict))
	infile = open(infn, 'r')
	for line in infile:
		line = line.rstrip()
		linedict = {}
		values = line.split("\t")
		if (values[2] == "exon"):
			#print line
			attributes = values[8].split(";")
			del attributes[-1]
			attributes = [x.strip() for x in attributes]
			for attribute in attributes:
				attr = attribute.strip().split(" ")
				linedict[attr[0]] = attr[1].strip("\"")
			txedge = linedict['transcript_id'] + "_" + values[3]
			txedgedict[txedge] = values[4]		
			#print txedge
			txedge = linedict['transcript_id'] + "_" + values[4]
			txedgedict[txedge] = values[3]		
			#print txedge
	infile.close()
	return txedgedict


def temp(line,lookup,txedgedict):
	"""
	Classify event types
	Identify inclusion/exclusion paths
	Generate "joinstrings"
	Add geneid using lookup table
	"""
	linedict = {}
	flanks1 = []
	flanks2 = []
	chains1 = []
	chains2 = []

	chr = line['seqname']
	astacode = line['structure']
	strand = line['strand']
	transcripts = line['transcript_id']
	astagene = line['gene_id']
	"""
	Get geneid from txids
	"""
	txie = transcripts.split(',')
	txlist1 = txie[0].split('/')
	txlist2 = txie[1].split('/')
	txlist = txlist1
	txlist.extend(txlist2)

	eventcode = astacode_to_english(astacode)
	flankstring = line['flanks']
	chainstring = line['splice_chain']
			
	joinstring = chain_to_exons(flankstring,chainstring,chr,strand,eventcode)

	print "FOO:", joinstring


	#null,22413976) 22415095(22414879^22414817-,22415052(

	digits = re.compile('\d+|null')
	#m = pattern.search(flanks[0])
	
	chains = chainstring.split(",");  # List of sites
	flanks = flankstring.split(",");  # The flanks (Always 2?)

	flanks1 = re.findall('\d+|null', flanks[0])
	flanks2 = re.findall('\d+|null', flanks[1])
	
	chains1 = re.findall('\d+', chains[0])
	chains2 = re.findall('\d+', chains[1])
	
	#print "**", flankstring, chainstring
	#print ",".join(flanks1)
	#print ",".join(flanks2)
	#print ",".join(chains1)
	#print ",".join(chains2)
	#print len(chains2)

	# For each event, get the start and end coordinates for the unit
	# this will be the connection to each edge in flanks.
	if (flanks1[0] == "null"):
		mykey = txlist1[0] + "_" + chains1[0]
	else:
		mykey = txlist1[0] + "_" + flanks1[0]
	left = txedgedict[mykey]
	#print "L:", left, mykey

	if (flanks2[0] == "null"):
		mykey = txlist2[0] + "_" + chains2[-1]
	else:
		mykey = txlist2[0] + "_" + flanks2[0]
	right = txedgedict[mykey]
	#print "R:", right, mykey
	
	eventstart = min(left,right)
	eventend = max(left,right)
	#print eventstart, eventend
	
	elements = [str(chr),str(eventstart),str(eventend),str(strand)]
	coreid = ":".join(elements)
	#print coreid
		

def asta_gff3(line,lookup):
	"""
	Classify event types
	Identify inclusion/exclusion paths
	Generate "joinstrings"
	Add geneid using lookup table
	"""
	linedict = {}

	digits = re.compile('\d+|null')
	#m = pattern.search(flanks[0])
	

	flankstring = line['flanks']
	chainstring = line['splice_chain']
	
	flankcoords1 = NULL
	
	
	chr = line['seqname']
	astacode = line['structure']
	strand = line['strand']
	transcripts = line['transcript_id']
	astagene = line['gene_id']
	"""
	Get geneid from txids
	"""
	txie = transcripts.split(',')
	txlist1 = txie[0].split('/')
	txlist2 = txie[1].split('/')
	txlist = txlist1
	txlist.extend(txlist2)
	gnames, gids = geneidsFromTxlist(txlist,lookup)
	
	eventcode = astacode_to_english(astacode)
			
	joinstring = chain_to_joins(flankstring,chainstring,chr,strand,eventcode)

	"""
	Excluding some events
	(Be conservative)
	"""
	
	if (astacode == "0,0") or (eventcode == "FELE"):
		pass
	else:
		'''
		Instantiate events
		'''
		chains = chainstring.split(",");  # List of sites
		flanks = flankstring.split(",");  # The flanks (Always 2?)
		print line
		# Get outer edges of event	
		print flanks[0], flanks[1]
		pattern = re.compile('\d+')
		m = pattern.search(flanks[0])
		lflank = m.group(0)
		m = pattern.search(flanks[1])
		rflank = m.group(0)
		
		print lflank, rflank
		if flanks[0] == "null":
			quit("NUll")
		if int(rflank) < int(lflank):
			quit("Order error")
		

		linedict = line
		linedict['joinstring'] = joinstring
		joins = joinstring.split(";")
		joinstring1 = joins[0]
		joinstring2 = joins[1]
		joins1 = joinstring1.split(',')
		joins2 = joinstring2.split(',')
		linedict['txlist1'] = txlist1
		linedict['txlist2'] = txlist2
		linedict['eventcode'] = eventcode
		linedict['gnames'] = gnames
		linedict['gids'] = gids



	return linedict
    
def main():
	
		
 	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 	# Intializing the reference
 	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 	# You need the gtf file, and the fasta file
 	#global lookup
 	#lookup = prep_ref_sub.prep_ref(gtffile,fastafile)
  	## Note that you now have a reference called ref.bam, and a lookup dict
 	#reffile = "ref.bam"
 
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 	# Intializing the reference
 	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 	# You need the gtf file, and the fasta file
 	#lookup = spanky_utils.prep_ref(gtffile,fastafile,output_dir)
  	## Note that you now have a reference called ref.bam, and a lookup dict
	#tmp_dir = output_dir + "/tmp/"
 	#reffile = tmp_dir + "/ref.bam"
	gtffile = "/users/davidsturgill/data/annotation/current.gtf"

	txdict = spanky_utils.gtf_to_attributes_dict(gtffile)
	lookup = txdict
	
	#for txid in txdict:
	#	print txid, txdict[txid]

	txedgedict = gtf_to_tx_edges_dict(gtffile)

	#for txedge in txedgedict:
	#	print txedge, txedgedict[txedge]

 	########################################
	########################################
	# Loading in data
	########################################
	########################################
 	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 	# Loading asta defs
 	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	astafile = "/users/davidsturgill/data/annotation/splicing_defs/astalavista_BDGP5.25.62_k2.out"

 	
	#print >> sys.stderr, "[%s] Loading asta defs %s" % (spanky_utils.timestamp(), astafile)
	astalines = GTFtoDict(astafile)

	# data look like this:	
	#{'seqname': 'chr2L', 'end': '9484', 'localization': 'CDS-3UTR_maxCDS', 'degree': '4', 'start': '8193', 'unknown': '.', 'flanks': '8193-,9484)', 'feature': 'as_event', 'gene_id': 'chr2L:7529-9484W', 'source': 'Undefined', 'score': '2.00000', 'sources': 'Undefined,Undefined', 'structure': '0,1^2-', 'splice_chain': ',8589^8668-', 'attribute': 'transcript_id "FBtr0300689,FBtr0300690"; gene_id "chr2L:7529-9484W"; flanks "8193-,9484)"; structure "0,1^2-"; splice_chain ",8589^8668-"; sources "Undefined,Undefined"; dimension "2_2"; degree "4"; localization "CDS-3UTR_maxCDS"; \n', 'transcript_id': 'FBtr0300689,FBtr0300690', 'dimension': '2_2', 'strand': '+'}
	

 	########################################
	########################################
	# Analyses
	########################################
	########################################

 	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 	# Build intron chains
 	# Get coordinates to measure
 	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 	
	print >> sys.stderr, "[%s] Iterating through events" % (spanky_utils.timestamp())
	
	# Try iterating over lines, calling sub many times
	eventdict = collections.defaultdict(lambda : collections.defaultdict(dict))	
	counter = 0
	excluded = 0
	for line in astalines:
		linedict = parse_asta_defs(line,lookup)
		#temp(line,lookup,txedgedict)
		if linedict:
			counter += 1		
			eventid = 'ASTA' + "%05d" % (counter,)
			eventdict[eventid] = linedict
			#print eventid, eventdict[eventid]['joins1'], eventdict[eventid]['joins2']
			chr = linedict['seqname']
			astacode = linedict['structure']
			strand = linedict['strand']
			structure = linedict['structure']

			joins1 = []
			for junc in eventdict[eventid]['joins1']:
   				coords = junc.split(":")
   				try:
   					joins1.extend(coords[1].split("_"))
   				except:
   					pass
			joins2 = []
			for junc in eventdict[eventid]['joins2']:
   				coords = junc.split(":")
   				try:
   					joins2.extend(coords[1].split("_"))
   				except:
   					pass
			#if (joins1 and joins2):
			if (eventdict[eventid]['structure'] == "0,1-2^"):
				coords1 = joins1
				coords2 = joins2

				tx1 = eventdict[eventid]['txlist1'][0]
				tx2 = eventdict[eventid]['txlist2'][0]
	
				edge = int(min(coords1)) - 1
				l1_key = tx1 + "_" + str(edge)
				edge = int(max(coords1)) + 1
				r1_key = tx1 + "_" + str(edge)
				#print l1_key
				#print r1_key
				l1 = txedgedict[l1_key]
				r1 = txedgedict[r1_key]
				#print l1, eventdict[eventid]['joins1'], r1
	
				edge = int(min(coords2)) - 1
				l2_key = tx2 + "_" + str(edge)
				edge = int(max(coords2)) + 1
				r2_key = tx2 + "_" + str(edge)
				l2 = txedgedict[l2_key]
				r2 = txedgedict[r2_key]
				#print l2, eventdict[eventid]['joins2'], r2
				
				# Build an exon chain
				featurestart = l1
				featureend = r1
				mycoords = [featurestart]
				mycoords.extend(coords2)
				mycoords.append(featureend)
				
				#print "mycoords", len(mycoords)
				#print mycoords
				#print mycoords[0]
				#print mycoords[1]
				#print mycoords[2]
				#print mycoords[3]
				exons = []
				exonelements = [str(chr), str(mycoords[0]), str(mycoords[1]), str(strand)]
				exons.append(":".join(exonelements)) 			
				exonelements = [str(chr), str(mycoords[2]), str(mycoords[3]), str(strand)]
				exons.append(":".join(exonelements)) 			
				exonelements = [str(chr), str(mycoords[4]), str(mycoords[5]), str(strand)]
				exons.append(":".join(exonelements)) 			

				#print "".join(exons)
				ID = "@".join(exons)
				#print len(exons)
				#print "ID:", ID

				line = []
				line.append(chr)
				line.append(structure)
				line.append("gene")
				line.append(featurestart)
				line.append(featureend)
				line.append(".")
				line.append(strand)
				line.append(".")
				featureline = "ID=" + ID + ";Parent=" + ID
				line.append(featureline)
				print "\t".join(line)

				subid = ".A"

				line = []
				line.append(chr)
				line.append(structure)
				line.append("mRNA")
				line.append(featurestart)
				line.append(featureend)
				line.append(".")
				line.append(strand)
				line.append(".")
				featureline = "ID=" + ID + subid + ";Parent=" + ID
				line.append(featureline)
				print "\t".join(line)

				subid = ".B"
				
				line = []
				line.append(chr)
				line.append(structure)
				line.append("mRNA")
				line.append(featurestart)
				line.append(featureend)
				line.append(".")
				line.append(strand)
				line.append(".")
				featureline = "ID=" + ID + subid + ";Parent=" + ID
				line.append(featureline)
				print "\t".join(line)

				line = []
				line.append(chr)
				line.append(structure)
				line.append("exon")
				line.append(featurestart)
				line.append(featureend)
				line.append(".")
				line.append(strand)
				line.append(".")
				featureline = "ID=" + ID + subid + ";Parent=" + ID
				line.append(featureline)
				print "\t".join(line)

				#21082714 ['chrX:21082827_21150208:-'] 21150213
				#21082714 ['chrX:21094675_21150208:-', 'chrX:21082827_21094606:-'] 21150213

				
		else:
			excluded += 1
	quit("Done")
	print "Processed", counter + excluded, "events"
	print "Excluded", excluded, "events"

	print >> sys.stderr, "[%s] Done processing events %s" % (spanky_utils.timestamp(), output_dir)
	
	summarize_events(eventdict,asta_summary_out)

	# A sample entry
	#{'attribute': 'transcript_id "FBtr0089158,FBtr0089162"; gene_id "chr4:89956-129430W"; flanks "null,92947-"; structure "1(2^,3(4^"; splice_chain "89956(90020^,91032(91193^";
	#sources "Undefined,Undefined"; dimension "2_6"; degree "4"; localization "5UTR_maxCDS"; \n', 'joins2': ['chr4:91194_92946:+'], 'joins1': ['chr4:90021_92946:+'], 
	#'sources': 'Undefined,Undefined', 'splice_chain': '89956(90020^,91032(91193^', 'seqname': 'chr4', 'end': '92947', 'eventcode': 'altFE', 'source': 'Undefined', 'unknown': '.',
	#'feature': 'as_event', 'start': '89956', 'score': '1.22222', 'strand': '+', 'gnames': 'FBgn0085432', 'degree': '4', 'localization': '5UTR_maxCDS', 'gene_id': 'chr4:89956-129430W', 
	#'transcript_id': 'FBtr0089158,FBtr0089162', 'flanks': 'null,92947-', 'structure': '1(2^,3(4^', 'gids': 'pan', 'txlist2': ['FBtr0089162'], 'txlist1': ['FBtr0089158', 'FBtr0089162'], 
	#'joinstring': 'chr4:90021_92946:+;chr4:91194_92946:+', 'dimension': '2_6'}


			
if __name__ == "__main__":
    sys.exit(main())





