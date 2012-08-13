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
	parser = argparse.ArgumentParser(description='Parse a bam file.')
	parser.add_argument('-i', help='bam file name', action="store", dest="i")
	parser.add_argument('-a', help='anchor size', action="store", type=int, dest="a", default=8)
	parser.add_argument('-jlist', help='junctionlist', action="store", dest="jlist")
	parser.add_argument('-jtab', help='junctiontab', action="store", dest="jtab")
	parser.add_argument('-o', help='outfile', action="store", dest="o", default="curated_juncs")
	parser.add_argument('-f', help='Fasta file', action="store", dest="f")
	parser.add_argument('-gff', help='GFF file', action="store", dest="gff")
	parser.add_argument('-gtf', help='GTF file', action="store", dest="gtf")
	parser.add_argument('-juncbed', help='BED file', action="store", dest="juncbed")
	parser.add_argument('-intronbed', help='BED file', action="store", dest="intronbed")
	
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
jtools

Performs operations on junctions, and converts between various formats

-----------------------------------------------------------------
'''

args = parse_options(desc)
bamfile = args.i
anchorsize = args.a
jlist = args.jlist
jtab = args.jtab
outfile = args.o
fastafile = args.f
juncbedfile = args.juncbed
intronbedfile = args.intronbed
gfffile = args.gff
gtffile = args.gtf

#~~~~~~~~~~~~~~~~~~~
# Prepare output directory
#~~~~~~~~~~~~~~~~~~~
output_dir = "."
#spanki_utils.prepare_basic_output_dir(output_dir)
spanki_utils.prepare_output_dir(output_dir)


#~~~~~~~~~~~~~~~~~~~
# Prepare output file
#~~~~~~~~~~~~~~~~~~~
# Prepare output file names
#juncs_out_name = output_dir + "/juncs.all"
#juncs_out = open(juncs_out_name, "w")

def main():
	
	# Functions for jtools:
	# Get sequence
	# Detect tandem acceptors NAGNAG
	# Annotate with genes
	# Jiggle
	# Bed to juncid
	# Guess frame
	# Find stops in intron + in frame
	# SVM recomputes
	# Splice site strength? ppt? 

	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# Load a fasta file
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	global f
	# Opening fasta filehandle
	print >> sys.stderr, "[%s] Opening fasta file" % (spanki_utils.timestamp())
	f = Fasta(fastafile)
	fastachr = set(sorted(f.keys()))
	#print fastachr

	########################################################
	### Parsing a juncbed file
	########################################################
	if (juncbedfile):

		print "juncid\toriginal_id\tdastring"
		print >> sys.stderr, "Loading", juncbedfile
		lines = csv.reader(open(juncbedfile, 'rb'), delimiter='\t')
		z = []
		for line in lines:
			pattern = re.compile('track')
			track = pattern.search(line[0])
			if not track:
				values = line
				blocksizes = values[10].split(",")
				blockstarts = values[11].split(",")
				chr = values[0]
				rangestart = int(values[1]) - 1
				rangeend = int(values[2])
				strand = values[5]
				id = values[3]
				intronstart = rangestart + int(blocksizes[0]) + 2
				intronend = rangeend - int(blocksizes[1]) 
				# Or..
				#intronend = rangestart + int(blocksizes[0]) + int(blockstarts[1])
			
				#chrXHet	800	1767	JUNC00000001	2	+	800	1767	255,0,0	2	20,63	0,904
			
				intronsize = intronend - intronstart;
			
				juncid = chr + ":" + str(intronstart) + "_" + str(intronend) + ":" + strand
				dastring = intron_sequence_single(juncid,f)
				z.append(str(dastring))
				print juncid, values[3], dastring
		
		print >> sys.stderr, "Distribution of detected motifs:\n",Counter(z)
		quit("Done")
	########################################################
	### Parsing a intronbed file
	########################################################
	#scaffold_12916	13833982	13834044	10
	#scaffold_12916	13838614	13838676	67
	#scaffold_12916	13839119	13839204	75

	if (intronbedfile):
		print "juncid\tid\tdastring"
		lines = csv.reader(open(intronbedfile, 'rb'), delimiter='\t')
		for line in lines:
			pattern = re.compile('track')
			track = pattern.search(line[0])
			values = line
			if not track:
				chr = values[0]
				intronstart = int(values[1]) + 1
				intronend = int(values[2]) - 1
				strand = "+"
				id = values[0]
			
				intronsize = intronend - intronstart;
			
				juncid = chr + ":" + str(intronstart) + "_" + str(intronend) + ":" + strand
				dastring = intron_sequence_single(juncid,f)
				
				print juncid, values[3], dastring
		
		
		quit("Done")
	########################################################

	########################################################
	### Converting from another format
	########################################################

	if gfffile:
	
		#reflist = tab_to_dict(gff)
		results = collections.defaultdict(lambda : collections.defaultdict(dict))
		gffdict = gff_to_dict(gfffile)
		for x in gffdict:
			#print x
			#print gffdict[x]
			if (gffdict[x]['feature_type'] == "exon_junction"):
				juncid = gffdict[x]['chr'] + ":" + str(int(gffdict[x]['start']) + 1) + "_" + str(int(gffdict[x]['end']) - 1) + ":" + gffdict[x]['strand']
			elif (gffdict[x]['feature_type'] == "intron"):
				juncid = gffdict[x]['chr'] + ":" + gffdict[x]['start'] + "_" + gffdict[x]['end'] + ":" + gffdict[x]['strand']
			dastring = intron_sequence_single(juncid,f)
			#print dastring
			results[x]['juncid'] = juncid
			results[x]['dastring'] = dastring
	
		print "ID\tjuncid\tdastring"
		for x in sorted(results.iterkeys()):
			print x, "\t", results[x]['juncid'], "\t", results[x]['dastring']
	
		quit()

	########################################################
	### Converting from another format
	########################################################

	if gtffile:
	
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# Intializing the reference
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# You need the gtf file, and the fasta file
		lookup = spanki_utils.prep_ref(gtffile,fastafile,output_dir)
		## Note that you now have a reference called ref.bam, and a lookup dict
		#tmp_dir = output_dir + "/tmp/"
		#reffile = tmp_dir + "/ref.bam"
		reffile = "tmp/ref.bam"
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# Load an annotation, flattened as bam
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		print >> sys.stderr, "[%s] Trying to load annotation as bam" % (spanki_utils.timestamp())
		reffh = pysam.Samfile( reffile, "rb" )
		edgedict, refjuncs = spanki_parse_utils.parseRefAsBam(reffh)
		reffh.close()
		print >> sys.stderr, "[%s] Done loading annotation as bam" % (spanki_utils.timestamp())
	
		for junc in refjuncs:
			print junc

	
		quit()
	### Below are functions that operate on a junction list
	########################################################


	if jlist:
	
		#~~~~~~~~~~~~~~~~~~~
		# Load reference junction list
		#~~~~~~~~~~~~~~~~~~~
		reflist = tab_to_dict(jlist)
	
		# Find the junctions in jlist that are not in jtab
	
		myjuncs = reflist.keys()
	
	
	
		
		print >> sys.stderr, len(myjuncs), "in junction list"
	
		updonor = 20
		downdonor = 2
		upacceptor = 2
		downacceptor = 20
		
	
		for x in myjuncs:
			print x
			j1 = Junctionid(x)
			j1.display()
			if j1.strand == "+":
				#print Seq(f[j1.chr][j1.donor-updonor:j1.donor], IUPAC.unambiguous_dna)
				tempseq = Seq(f[j1.chr][j1.donor-updonor:j1.donor], IUPAC.unambiguous_dna)
				#print "***", tempseq.translate()
				
				#print Seq(f[j1.chr][j1.donor:j1.donor + downdonor], IUPAC.unambiguous_dna)
				#print Seq(f[j1.chr][j1.acceptor-upacceptor:j1.acceptor], IUPAC.unambiguous_dna)
				#print Seq(f[j1.chr][j1.acceptor:j1.acceptor + downacceptor], IUPAC.unambiguous_dna)
				nagstring = find_nag(Seq(f[j1.chr][j1.acceptor:j1.acceptor + downacceptor], IUPAC.unambiguous_dna))
				print nagstring
			elif j1.strand == "-":
				pass
				#print Seq(f[j1.chr][j1.donor:j1.donor + updonor], IUPAC.unambiguous_dna).reverse_complement()
				#print Seq(f[j1.chr][j1.donor - downdonor:j1.donor], IUPAC.unambiguous_dna).reverse_complement()
				#print Seq(f[j1.chr][j1.acceptor:j1.acceptor + upacceptor], IUPAC.unambiguous_dna).reverse_complement()
				#print Seq(f[j1.chr][j1.acceptor-downacceptor:j1.acceptor], IUPAC.unambiguous_dna).reverse_complement()
			else:
				quit("Don't recognize strand")
	
				#fiveprimeflank = fiveprimeflank.reverse_complement()
		quit("Done")
	
	
	quit()
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# Older code that's not used yet
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
	
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# IRT
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	bamfh = pysam.Samfile( bamfile, "rb" )
	#for alignedread in samfile:
	# Need some kind of iterator to getread length from first alignment in sam
	print >> sys.stderr, "[%s] Getting intron read-though (IRT), may take awhile" % (spanki_utils.timestamp())
	IRT = intron_readthrough(myjuncs,bamfh)
	bamfh.close()
	print >> sys.stderr, "[%s] Done getting IRT" % (spanki_utils.timestamp())



	#for edgeid in covbyedge.keys():
	#	print edgeid, covbyedge[edgeid]
		

	# These are the fields you end up with after merging:
	#juncid	geneassign	cov	lirt	rirt	irt	dncov	ancov	numsamps
	#chr2L:22427471_22427525:- 	none 	2 	57 	28 	85 	0 	0 	1
	#chr2R:5702257_5702656:+ 	FBgn0040092 	13 	0 	0 	0 	0 	0 	2
	#chr2L:11436293_11436415:- 	FBgn0261648 	23 	0 	0 	0 	0 	0 	2
	#chr2R:9334834_9336812:- 	FBgn0013765 	6 	0 	0 	0 	0 	0 	2

	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# Now compile the results
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# First show how you can get in hte myjuncs list
	print >> sys.stderr, "Printing results table"
	print >> juncs_out, "juncid\tgeneassign\tannostatus\tintron_size\tgmcode\tregcode\tcov\tlirt\trirt\tirt\tdncov\tancov"
	
	for juncid in sorted(keys2):
		try:
			results = [juncid, jdict[juncid]['geneassign'], jdict[juncid]['annostatus'], jdict[juncid]['intron_size'], jdict[juncid]['gmcode'], jdict[juncid]['regcode'], jdict[juncid]['cov'], jdict[juncid]['lirt'], jdict[juncid]['rirt'], jdict[juncid]['irt'], jdict[juncid]['dncov'], jdict[juncid]['ancov']]
			print >> juncs_out, ('\t'.join(map(str,results)))
		except KeyError:
			#myjuncs.append(juncid)	
			j1 = Junctionid(juncid)
			donid = j1.donid
			accid = j1.accid
			if covbyedge[donid]: dncov = covbyedge[donid]
			else: dncov = 0
			if covbyedge[accid]: ancov = covbyedge[accid]
			else: ancov = 0
			results = [juncid, reflist[juncid]['geneassign'],  reflist[juncid]['annostatus'],  reflist[juncid]['intron_size'],  reflist[juncid]['gmcode'],  reflist[juncid]['regcode'], 0, IRT[juncid]['lirt'], IRT[juncid]['rirt'], IRT[juncid]['irt'], dncov, ancov]
			#print(results, sep='\t')
			print >> juncs_out, ('\t'.join(map(str,results)))

	quit("done")

 	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 	# Parse the read alignments
 	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# Parse the bam file
	## Get a table of junctions, table of donors etc.
	bamfh = pysam.Samfile( bamfile, "rb" )
	#JTAB,UNFILT_JTAB,STAB,NEWDTAB,MMES = parse_aligns_detailed(bamfh)
	JTAB,UNFILT_JTAB = quickcov(bamfh,anchorsize)
	bamfh.close()
	myjuncs = JTAB.keys()
	myjuncs.sort()
	
 	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 	# Print junction list to the output directory
 	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 	print "juncid\tunfilt_cov\tcov"
	for juncid in myjuncs:
		print juncid, UNFILT_JTAB[juncid], JTAB[juncid]

	
	

if __name__ == "__main__":
    sys.exit(main())

