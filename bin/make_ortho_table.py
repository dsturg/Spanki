#!/usr/bin/env python
# encoding: utf-8
"""
jtools

Performs various operations on junctions, and does
format conversions
"""
from __future__ import division 

import re
import sys
import argparse
import pysam
import csv
import math
import os
import collections
import operator

# For counting dastring types 
from collections import Counter

from pyfasta import Fasta
from datetime import datetime, date

# Biopython, for revcomp
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

# Custom modules to import:
import spanki.spanki_parse_utils as spanki_parse_utils
import spanki.spanki_utils as spanki_utils

class Junctionid:
	"""
	Base class for a junction, from juncid
	Expects junction id
	"""
	def __init__(self, juncid):
		chr = juncid.split(':')[0]
		coords = juncid.split(':')[1]
		strand = juncid.split(':')[2]
		start = int(coords.split('_')[0]) - 1
		end = int(coords.split('_')[1])
		self.chr = chr
		self.start = start
		self.end = end
		self.strand = strand.strip()
		self.intronsize = end - start
		#self.accid = str(coords.split('_')[0])
		if self.strand == "+":
			self.donid = chr + ":" + str(start)
			self.donor = start
			self.accid = chr + ":" + str(end)
			self.donor = start
			self.acceptor = end
		elif self.strand == "-":
			self.accid = chr + ":" + str(start)
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
			

class JuncBasic:
	"""
	Base class for a junction
	Expects SAM input
	"""
	def __init__(self, chr, start, cigar, readseq, tags):
		self.chr = chr
		self.start = start
		self.cigar = cigar
		self.readseq = readseq
		## anchorsize is the genomic coordinate space (accounting for indels)
		## The aligned portion may be longer or shorter after indels
		left_anchorsize = 0
		left_read_anchor = ""
		left_genome_aligned = ""
		right_anchorsize = 0
		right_read_anchor = ""
		right_genome_aligned = ""
		gapsize = 0
		gap_seen = 0
		readpos = 0
		genomepos = self.start
		for i in self.cigar:
			if i[0] == 0:
				if gap_seen < 1:
					left_anchorsize += i[1]
					left_read_anchor += self.readseq[readpos:readpos + i[1]]
					left_genome_aligned += "M" * (i[1] + 1) 
				else:
					right_anchorsize += i[1]
					right_read_anchor += self.readseq[readpos:readpos + i[1]]
					right_genome_aligned += "M" * (i[1] + 1) 
				readpos += i[1]
				genomepos += i[1]
			elif i[0] == 1:
				'''
				If insertion, add to readpos
				'''
				readpos += i[1]
			elif i[0] == 2:
				'''
				If deletion, add to anchor size
				'''
				if gap_seen < 1:
					left_anchorsize += i[1]
				else:
					right_anchorsize += i[1]
				genomepos += i[1]
			elif i[0] == 3:
				gap_seen += 1
				gapsize += i[1]
			else:
				quit("Don't recognize cigar code")

		junc_left_side = self.start + left_anchorsize + 1
		junc_right_side = self.start + left_anchorsize + gapsize
		self.junc_left_side = junc_left_side
		self.junc_right_side = junc_right_side
		
		self.left_genome_aligned = left_genome_aligned
		self.right_genome_aligned = right_genome_aligned
		self.left_read_anchor = left_read_anchor
		self.right_read_anchor = right_read_anchor
		self.gapsize = gapsize

	
		""" 
		Try to get strand from SAM file first
		If it can't do it, use donor/acceptor motifs
		If still no, use *
		"""
		keys = []
		values = []
		for tag in tags:
			keys.append(tag[0])
			values.append(tag[1])
		tagdict = dict(zip(keys, values))
		try:
			strand = tagdict['XS']
			self.strand = strand
		except:
			self.strand = "*"
		
		self.juncid = self.chr + ":" + str(junc_left_side) + "_" + str(junc_right_side) + ":" + self.strand

def subdivideCigar(cigar):
	# Get gap indices
	gaps = []
	counter = 0
	for i in cigar:
		if i[0] == 3: gaps.append(counter)
		counter += 1	
	if gaps:
		"""
		If there are no gaps(juncs)
		return empty set
		"""
		borders = [0]
		for i in range(len(gaps[:-1])):
			borders.append(gaps[i+1])
			borders.append(gaps[i] + 1)
		borders.append(len(cigar))
		
		# Define subcigar ranges
		subcigars = []
		myjuncs = zip(borders[0::2], borders[1::2])
		for junc in myjuncs:
			subcigars.append(cigar[junc[0]:junc[1]])
	else:
		subcigars = []
	return subcigars


def quickcov(samfile,anchorsize):
	"""
	Takes a sam file as input
	Quickly gets coverage
	"""
	JTAB = collections.defaultdict(int)
	UNFILT_JTAB = collections.defaultdict(int)
	zero_anchor_warnings = 0
	for alignedread in samfile:
		if (len(alignedread.cigar) > 1):
			# Note that alignedread.is_reverse does not work right to get strand
			#strand = alignedread.tags[1][1] 
			mytid = alignedread.tid
			chr = samfile.getrname(mytid)
			start = alignedread.pos
			offset = alignedread.pos
			cigar = alignedread.cigar
			readseq = alignedread.query

			subcigars = subdivideCigar(cigar)

			for cigar in subcigars:
				j1 = JuncBasic(chr,start,cigar,readseq,alignedread.tags)
				'''
				For multi-gap cigars, you have to add to the start
				for the next one
				'''
				start += int(len(j1.left_genome_aligned))
				start += j1.gapsize

				juncid = j1.juncid

				
				'''
				Apply a restriction on anchor size
				to count the junction
				'''
				min_anchor = min(len(j1.left_read_anchor),len(j1.right_read_anchor))

				'''
				Also count unfiltered (no anchor cutoff)
				'''

				if min_anchor < 1:
					'''
					Warn when anchor size is zero
					For example, cigar 76M54354N
					These are errors, and shouldn't count toward coverage
					'''
					zero_anchor_warnings += 1
				else: 
					UNFILT_JTAB[juncid] += 1; # Total junction spanning reads


				if min_anchor >= anchorsize:
				
					JTAB[juncid] += 1; # Total junction spanning reads
	
				else:
					pass

	if zero_anchor_warnings > 0:
		print "[ WARNING ] ", zero_anchor_warnings, "alignments with zero anchor coverage excluded"
	
	return JTAB,UNFILT_JTAB

def tab_to_dict(tabfile):
	"""
	Generic make a dict from a table
	Assumes first column has key
	and there are column headers
	"""
	mytab = {}
	lines = csv.reader(open(tabfile, 'rb'), delimiter='\t')
	linecount = 0
	for line in lines:
		if (linecount < 1):
			"""
			First line is column header - use as keys
			"""
			keys = line
		else: 
			values = line
			linedict = dict(zip(keys, values))
			id = str(values[0])
			mytab[id] = linedict 
			#print "adding to ",linedict['juncid']
			#print linedict
		linecount += 1
	return mytab

def intron_readthrough(myjuncs,samfile):
	IRT = collections.defaultdict(lambda : collections.defaultdict(dict))
	overhang = 8
	'''
	Get read length from first aligment
	'''
	for alignedread in samfile:
		readlength = alignedread.rlen
		break
	for juncid in myjuncs:

		j1 = Junctionid(juncid)

		# Check left side
		lirt = 0
		rangestart = j1.start - readlength - overhang
		rangeend = j1.start - overhang
		if (rangestart < 0): rangestart = 0
		# Need this for cases where a junciton is near 
		# the end of a chromosome
		
		
		# Get all reads in the region where intron read-thru reads
		# could reside.
		# Count all reads that start in the correct range and have no gaps.
		reads = samfile.fetch( j1.chr, rangestart, rangeend )		
		for read in reads:
			if (rangeend > read.pos > rangestart): 
				cigar = read.cigar
				gaps = 0
				if len(cigar) > 1:
					for i in cigar:
						if i[0] == 3: gaps += 1
				if gaps < 1: 
					lirt += 1
		# Check right side
		rirt = 0
		rangestart = j1.end - readlength + overhang
		rangeend = j1.end - overhang
		if (rangestart < 0): rangestart = 0
		reads = samfile.fetch( j1.chr, rangestart, rangeend )
		for read in reads:
			if (rangeend > read.pos > rangestart): 
				cigar = read.cigar
				gaps = 0
				if len(cigar) > 1:
					for i in cigar:
						if i[0] == 3: gaps += 1
				if gaps < 1: 
					rirt += 1
		IRT[juncid]['lirt'] = lirt
		IRT[juncid]['rirt'] = rirt
		IRT[juncid]['irt'] = lirt + rirt
	return IRT

def donor_acceptor_transition(NEWDTAB,myjuncs):
	"""
	Gets donor-acceptor transition probabilities 
	using junction coverage
	"""
	DATRANS = collections.defaultdict(lambda : collections.defaultdict(dict))
	covbyedge = collections.defaultdict(lambda : collections.defaultdict(dict))
	joinsbyedge = collections.defaultdict(lambda : collections.defaultdict(dict))
	Ks = NEWDTAB.keys()
	Ks.sort()
	'''
	First calculate coverage by edge (donor or acceptor)
	Remember to check if each join has been filtered out
	'''
	for x in Ks:
		#dcov = sum(NEWDTAB[x].values())
		#covbyedge[x] = dcov
		covbyedge[x] = 0
		joinsbyedge[x] = 0
		for y in NEWDTAB[x].keys():
			'''
			Get juncid
			'''
			juncid = str(x) + "_" + str(y)
			if (juncid + ":+") in myjuncs:
				juncid = juncid + ":+"
			elif (juncid + ":-") in myjuncs:
				juncid = juncid + ":-"
			else:
				pass
			''' 
			Check to see if it was filtered
			'''
			if juncid in myjuncs:
				covbyedge[x] += NEWDTAB[x][y]
				joinsbyedge[x] += 1
				#jcov = NEWDTAB[x][y]
				#print "Looking at ", x, "to", y, "with jcov", jcov
				if covbyedge[y]:
					covbyedge[y] += NEWDTAB[x][y]
					joinsbyedge[y] += 1
				else:
					covbyedge[y] = NEWDTAB[x][y]
					joinsbyedge[y] = 1
	'''
	Now instantiate hash of neighbor coverages keyed by juncid
	Also do DATRANS
	'''
	NEIGHBORCOV = collections.defaultdict(lambda : collections.defaultdict(dict))
	for x in Ks:
		for y in NEWDTAB[x].keys():
			juncid = str(x) + "_" + str(y)
			if (juncid + ":+") in myjuncs:
				juncid = juncid + ":+"
			elif (juncid + ":-") in myjuncs:
				juncid = juncid + ":-"
			else:
				pass

			if (juncid in myjuncs):
				#print "Looking at neighbor coverage"
				jcov = NEWDTAB[x][y]
				#print "covbyedge", covbyedge[x]
				try:
					trans = jcov/covbyedge[x]
				except ZeroDivisionError:
					trans = 0
				DATRANS[juncid]['datrans'] = "%.3f" % trans
				NEIGHBORCOV[juncid]['dncov'] = covbyedge[x] - jcov
				NEIGHBORCOV[juncid]['ancov'] = covbyedge[y] - jcov
				NEIGHBORCOV[juncid]['dnjoins'] = joinsbyedge[x]
				NEIGHBORCOV[juncid]['anjoins'] = joinsbyedge[y]
	return DATRANS, NEIGHBORCOV, joinsbyedge

def find_nag(myseq):
	#IRT = collections.defaultdict(lambda : collections.defaultdict(dict))
	'''
	Find nag acceptors
	'''
	nagstring = ""
	print myseq
	nagstarts = []
	nagseqs = []

	p = re.compile('[ACGT]AG')
	#nags = p.findall(str(myseq))
	for nag in p.finditer(str(myseq)):
		print nag.start(), nag.group()
		nagstarts.append(str(nag.start()))
		nagseqs.append(nag.group())
		#nagstring = nagstring + "," + nag

	#p = re.compile("[a-z]")
	#for m in p.finditer('a1b2c3d4'):
    #print m.start(), m.group()
	if (nagseqs):
		nagstring = ','.join(nagseqs) + ";" + ','.join(nagstarts)
	else:
		nagstring = "-"
	return nagstring

################################
################################
################################
################################
################################
# GFF
#chr2L	FlyBase	exon_junction	3473704	3473770	.	+	.	ID=FBsf0000139598;Name=Dmel:r5:2L:3473704:3473770:+;exonA=2L:3473641..3473704:+;exonB=2L:3473770..3473835:+
#chr2L	FlyBase	intron	3473705	3473769	.	+	.	ID=intron_FBgn0051956:3_FBgn0051956:4;Name=pgant4-in;Parent=FBtr0077523,FBtr0114524;parent_type=mRNA
# GTF
#chr2L	protein_coding	exon	8193	9484	.	+	.	 gene_id "FBgn0031208"; transcript_id "FBtr0300689"; exon_number "2"; gene_name "CG11023"; transcript_name "CG11023-RB";
# BED
#chr4_group3	1015441	1025375	JUNC00024565	2	-	1015441	1025375	16711680	2	46,46,	0,9888,
#chr4_group3	1015827	1024519	JUNC00024566	7	-	1015827	1024519	16711680	2	60,67,	0,8625,

def gff_to_dict(infn):
	'''
	Returns a dict of the attributes of a GFF
	'''
	attrdict = collections.defaultdict(lambda : collections.defaultdict(dict))
	infile = open(infn, 'r')
	for line in infile:
		line = line.rstrip()
		linedict = {}
		values = line.split("\t")
		if len(values) > 1:
			if (values[2] == "exon_junction" or "intron" ):
				#print line
				attributes = values[8].split(";")
				del attributes[-1]
				attributes = [x.strip() for x in attributes]
				for attribute in attributes:
					attr = attribute.strip().split("=")
					print linedict
					linedict[attr[0]] = attr[1].strip("\"")
					#print "Keying on", attr[0]
				myfeatureid = linedict['ID']
				try:
					attrdict[linedict['ID']]['Name'] = linedict['Name']
				except:
					attrdict[linedict['ID']]['Name'] = ""
				attrdict[myfeatureid]['chr'] = values[0]
				attrdict[myfeatureid]['feature_type'] = values[2]
				attrdict[myfeatureid]['start'] = values[3]
				attrdict[myfeatureid]['end'] = values[4]
				attrdict[myfeatureid]['strand'] = values[6]
	infile.close()
	return attrdict

def intron_sequence(myjuncs,f):
 	"""
 	Returns the intron sequence and flanks for each join
 	"""
 	INTSEQ = collections.defaultdict(lambda : collections.defaultdict(dict))
	for juncid in myjuncs:
		j1 = Junctionid(juncid)
		if j1.strand == '+':
			fiveprimeflank = Seq(f[j1.chr][j1.start-10:j1.start], IUPAC.unambiguous_dna)
			threeprimeflank = Seq(f[j1.chr][j1.end:j1.end+10], IUPAC.unambiguous_dna)
			donormotif = Seq(f[j1.chr][j1.start:j1.start+2], IUPAC.unambiguous_dna)
			acceptormotif = Seq(f[j1.chr][j1.end-2:j1.end], IUPAC.unambiguous_dna)
			acceptormotif = acceptormotif.upper()
			donormotif = donormotif.upper()
			dastring = donormotif + '..' + acceptormotif
		else:
			fiveprimeflank = Seq(f[j1.chr][j1.end:j1.end+10], IUPAC.unambiguous_dna)
			threeprimeflank = Seq(f[j1.chr][j1.start-10:j1.start], IUPAC.unambiguous_dna)
			fiveprimeflank = fiveprimeflank.reverse_complement()
			threeprimeflank = threeprimeflank.reverse_complement()
			acceptormotif = Seq(f[j1.chr][j1.start:j1.start+2], IUPAC.unambiguous_dna)
			donormotif = Seq(f[j1.chr][j1.end-2:j1.end], IUPAC.unambiguous_dna)
			acceptormotif = acceptormotif.upper()
			donormotif = donormotif.upper()
			dastring = donormotif.reverse_complement() + '..' + acceptormotif.reverse_complement()
		INTSEQ[juncid]['dinucleotide'] = dastring
		INTSEQ[juncid]['flank5'] = fiveprimeflank
		INTSEQ[juncid]['flank3'] = threeprimeflank
	return INTSEQ

def intron_sequence_single(juncid,f):
 	"""
 	Returns the intron sequence and flanks for a join
 	"""
	j1 = Junctionid(juncid)
	if j1.strand == '+':
		fiveprimeflank = Seq(f[j1.chr][j1.start-10:j1.start], IUPAC.unambiguous_dna)
		threeprimeflank = Seq(f[j1.chr][j1.end:j1.end+10], IUPAC.unambiguous_dna)
		donormotif = Seq(f[j1.chr][j1.start:j1.start+2], IUPAC.unambiguous_dna)
		acceptormotif = Seq(f[j1.chr][j1.end-2:j1.end], IUPAC.unambiguous_dna)
		acceptormotif = acceptormotif.upper()
		donormotif = donormotif.upper()
		dastring = donormotif + '..' + acceptormotif
	else:
		fiveprimeflank = Seq(f[j1.chr][j1.end:j1.end+10], IUPAC.unambiguous_dna)
		threeprimeflank = Seq(f[j1.chr][j1.start-10:j1.start], IUPAC.unambiguous_dna)
		fiveprimeflank = fiveprimeflank.reverse_complement()
		threeprimeflank = threeprimeflank.reverse_complement()
		acceptormotif = Seq(f[j1.chr][j1.start:j1.start+2], IUPAC.unambiguous_dna)
		donormotif = Seq(f[j1.chr][j1.end-2:j1.end], IUPAC.unambiguous_dna)
		acceptormotif = acceptormotif.upper()
		donormotif = donormotif.upper()
		dastring = donormotif.reverse_complement() + '..' + acceptormotif.reverse_complement()
	#INTSEQ[juncid]['dinucleotide'] = dastring
	#INTSEQ[juncid]['flank5'] = fiveprimeflank
	#INTSEQ[juncid]['flank3'] = threeprimeflank
	return dastring

def print_dict_sorted_ordered(mydict,fout,sortfield,fieldorder,keyname):
	mykeys = fieldorder
	mytup = []
	print >> fout, keyname, "\t", '\t'.join(fieldorder)
	# Get tuples of event and sort field
	for x in mydict.keys():
		mytup.append([x,mydict[x][sortfield]])
	# Get list of keys sorted by sort field
	#sorted_mytup = sorted(mytup.iteritems(), key=operator.itemgetter(1))
	#sorted(student_tuples, key=itemgetter(2))
	sorted_keys = sorted(mytup, key=operator.itemgetter(1), reverse=False)
	mykeys = tuple(x[0] for x in sorted_keys)
	for x in mykeys:
		vals = []
		for field in fieldorder:
			vals.append(mydict[x][field])
		print >> fout, x, '\t', '\t'.join(map(str,vals))

def print_dict_sorted_ordered_stdout(mydict,sortfield,fieldorder,keyname):
	mykeys = fieldorder
	mytup = []
	print keyname, "\t", '\t'.join(fieldorder)
	# Get tuples of event and sort field
	for x in mydict.keys():
		mytup.append([x,mydict[x][sortfield]])
	# Get list of keys sorted by sort field
	#sorted_mytup = sorted(mytup.iteritems(), key=operator.itemgetter(1))
	#sorted(student_tuples, key=itemgetter(2))
	sorted_keys = sorted(mytup, key=operator.itemgetter(1), reverse=False)
	mykeys = tuple(x[0] for x in sorted_keys)
	for x in mykeys:
		vals = []
		for field in fieldorder:
			vals.append(mydict[x][field])
		print x, '\t', '\t'.join(map(str,vals))

################################
################################
################################
################################
################################
class MyParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        #sys.stderr.write(desc)
        self.print_help()
        sys.exit(2)

def parse_options(desc):
	parser=MyParser(description=desc, formatter_class=argparse.RawTextHelpFormatter)
	parser = argparse.ArgumentParser(description='Argument parser.')
	parser.add_argument('-i', help='input file', action="store", dest="i")
	
	args = parser.parse_args()
	# See here for argparse instructions:
	# http://docs.python.org/dev/library/argparse.html
	if len(sys.argv)==1:
		parser.print_help()
		sys.exit(1)
	return args

# Initialize parameters
desc = '''
-----------------------------------------------------------------
make_ortho_table.py

Makes a table from OrthoDB raw input
-----------------------------------------------------------------
'''

args = parse_options(desc)
inputfile = args.i

def main():
	'''
	Table looks like this:
	Drosophila	EOG5001BW	FBgn0130665	DGRIM	B4JSY1
	Drosophila	EOG5001BW	FBgn0133075	DMOJA	B4KE24
	Drosophila	EOG5001BW	FBgn0197452	DVIRI	B4M4V0
	Drosophila	EOG5001BW	FBgn0214180	DWILL	B4N969
	Drosophila	EOG5001BW	FBgn0150296	DPERS	B4GM15
	Drosophila	EOG5001BW	FBgn0077652	DPSEU	Q293J4
	Drosophila	EOG5001BW	FBgn0095545	DANAN	B3M2X7
	Drosophila	EOG5001BW	FBgn0107258	DEREC	B3P2Z8
	Drosophila	EOG5001BW	FBgn0242113	DYAKU	B4PQ99
	Drosophila	EOG5001BW	FBgn0013812	DMELA	Q9VDG0
	Drosophila	EOG5001BW	FBgn0170034	DSECH	B4IKM6
	Drosophila	EOG5001BW	FBgn0186760	DSIMU	B4NTE4
	Drosophila	EOG5001BW	FBgn0191501	DSIMU	B4R1N9
	Drosophila	EOG5001BX	FBgn0117540	DGRIM	B4JD53
	Drosophila	EOG5001BX	FBgn0139188	DMOJA	B4KFN4
	Drosophila	EOG5001BX	FBgn0203400	DVIRI	B4MDB2
	'''

	spporder = ['DMELA','DSIMU','DSECH','DYAKU','DEREC','DANAN','DPSEU','DPERS','DWILL','DMOJA','DVIRI','DGRIM']	


	'''
	The table I want
	
	#cluster_id	classification	dana	dere	dgri	dmel	dmoj	dper	dpse	dsec	dsim	dvir	dwil	dyak
	1	11111nn11111	dana_GLEANR_15484	dere_GLEANR_9302	dgri_GLEANR_10277	FBgn0031344	dmoj_GLEANR_1887	dper_GLEANR_1559,dper_GLEANR_1561	dpse_GLEANR_1674,dpse_GLEANR_1675	dsec_GLEANR_17788	dsim_GLEANR_6648	dvir_GLEANR_1064	dwil_GLEANR_15229	dyak_GLEANR_1708
	2	111111111111	dana_GLEANR_4448	dere_GLEANR_2763	dgri_GLEANR_2025	FBgn0030627	dmoj_GLEANR_6756	dper_GLEANR_17551	dpse_GLEANR_11073	dsec_GLEANR_25	dsim_GLEANR_17352	dvir_GLEANR_17450	dwil_GLEANR_16581	dyak_GLEANR_18528
	3	111111111110	dana_GLEANR_19778	dere_GLEANR_1506	dgri_GLEANR_747	FBgn0025865	dmoj_GLEANR_10244	dper_GLEANR_12820	dpse_GLEANR_6509	dsec_GLEANR_15893	dsim_GLEANR_3799	dvir_GLEANR_10103	dwil_GLEANR_13818	.
	4	11111n111111	dana_GLEANR_18417	dere_GLEANR_2219	dgri_GLEANR_1418	FBgn0037755	dmoj_GLEANR_9990	dper_GLEANR_5577,dper_GLEANR_5578	dpse_GLEANR_5476	dsec_GLEANR_9112	dsim_GLEANR_4532	dvir_GLEANR_9794	dwil_GLEANR_11359	dyak_GLEANR_8396
	5	111111111111	dana_GLEANR_7595	dere_GLEANR_2092	dgri_GLEANR_3125	FBgn0037958	dmoj_GLEANR_7867	dper_GLEANR_12529	dpse_GLEANR_6219	dsec_GLEANR_8979	dsim_GLEANR_4400	dvir_GLEANR_8751	dwil_G
	'''

	results = collections.defaultdict(lambda : collections.defaultdict(dict))

	print >> sys.stderr, "Loading", inputfile
	lines = csv.reader(open(inputfile, 'rb'), delimiter='\t')
	for line in lines:
		pattern = re.compile('#')
		track = pattern.search(line[0])
		if not track:
			values = line
			ogroup = values[1]
			fbgn = values[2]
			spp = values[3]
			#results[ogroup][spp] = fbgn
			if results[ogroup][spp]:
				results[ogroup][spp].append(fbgn)
			else:
				results[ogroup][spp] = []
				results[ogroup][spp].append(fbgn)
			#if spp == "DMELA":
			#	print "\t", ogroup, fbgn, spp

	'''
	Make a prettier table
	'''
	compiledtable = collections.defaultdict(lambda : collections.defaultdict(dict))
	print >> sys.stderr, "Processing table"
	for ogroup in results.keys():
		'''
		Make orthology string
		'''
		ortholist = ""
		orthostring = ""
		temp = []
		for species in spporder:
			if len(results[ogroup][species]) < 1:
				orthostring += "0"
			elif len(results[ogroup][species]) > 1:
				orthostring += "N"
			else:
				orthostring += "1"
			temp.append("%s" % len(results[ogroup][species]))
			if results[ogroup][species]:
				compiledtable[ogroup][species] = ",".join(results[ogroup][species])
			else:
				compiledtable[ogroup][species] = "-"
		#temp2 = map(chr, temp)
		compiledtable[ogroup]['ortholist'] = ",".join(temp)
		compiledtable[ogroup]['orthostring'] = orthostring
		#print temp
		#print results[ogroup]
		#for species in results[ogroup].keys
	#quit("Done")

	'''
	Output table
	'''
	#for ogroup in compiledtable.keys():
	#	print compiledtable[ogroup]
	#quit("Done")
	
	myresults = compiledtable
	fieldorder = ['orthostring','ortholist','DMELA','DSIMU','DSECH','DYAKU','DEREC','DANAN','DPSEU','DPERS','DWILL','DMOJA','DVIRI','DGRIM']	
	sortfield = 'ogroup'
	print_dict_sorted_ordered_stdout(myresults,sortfield,fieldorder,'ogroup')



	#for ogroup in results.keys():
	#	print results[ogroup]

	
	

if __name__ == "__main__":
    sys.exit(main())

