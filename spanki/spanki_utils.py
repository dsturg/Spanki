#!/usr/bin/env python
# encoding: utf-8

import collections
import subprocess
import sys
import os
import re
from datetime import datetime, date

def timestamp():
    currenttime = datetime.now()
    return currenttime.strftime("%b %d %H:%M:%S")

def gtf_to_attributes_dict(infn):
	'''
	Returns a dict of the attributes line of a GTF
	'''
	attrdict = collections.defaultdict(lambda : collections.defaultdict(dict))
	infile = open(infn, 'r')
	for line in infile:
		line = line.rstrip()
		linedict = {}
		values = line.split("\t")
		if len(values) > 1:
			if (values[2] == "exon"):
				#print "---------------\n"
				#print "line is\n";
				#print line
				#print "---------------\n"				
				attributes = values[8].split("; ")
				# Note there are gene names in Arabidopsis with a semicolon: "PIP1;2"
				# This will cause a parsing error here.
				# Note I now split on 'semicolon + space" - be on the lookout for exceptions
				# such as gene names with spaces
				del attributes[-1]
				attributes = [x.strip() for x in attributes]
				for attribute in attributes:
					attr = attribute.strip().split(" ") # Split each item in attibutes line by space
					try:
						linedict[attr[0]] = attr[1].strip("\"") # Remove quotes from field data, key by first string in attr
					except:
						print line # Print the offending line
						quit("GTF parsing error")
				try:
					attrdict[linedict['transcript_id']]['gene_id'] = linedict['gene_id']
				except:
					attrdict[linedict['transcript_id']]['gene_id'] = ""
				try:
					attrdict[linedict['transcript_id']]['gene_name'] = linedict['gene_name']
				except:
					attrdict[linedict['transcript_id']]['gene_name'] = "None"
	infile.close()
	return attrdict

# Example GTF line
#chr3R	protein_coding	exon	380	1913	.	+	.	 gene_id "FBgn0037213"; transcript_id "FBtr0078962"; exon_number "1"; gene_name "CG12581"; transcript_name "CG12581-RA";

def prep_ref(gtffile,fastafile,output_dir):
	'''
	From a gtf and fasta, creates a BAM representation
	Requries the gtf_to_sam function in Cufflinks (Trapnell et. al)
	http://cufflinks.cbcb.umd.edu
	'''
	tmp_dir = output_dir + "/tmp/"
	print >> sys.stderr, "[%s] Making transcript to attribute lookup" % (timestamp())
	txdict = gtf_to_attributes_dict(gtffile)
	print >> sys.stderr, "[%s] Convert GTF reference to SAM" % (timestamp())
	try:
		subprocess.call(["gtf_to_sam", gtffile, tmp_dir + "/ref.sam"])
	except:
		print "Error:  Can't call 'gtf_to_sam.'"  
		print "Please check that Cufflinks is installed and in your path."
		quit()		
	'''
	Check to see that fastafile is already indexed before trying to index it
	'''
	try:
		with open(fastafile + '.fai'): pass
	except IOError:
		try:
			subprocess.call(["samtools", "faidx", fastafile])
		except:
			print "Error:  Can't call samtools"
			print "Please check that samtools is installed"
			quit()

	fastidx = fastafile + ".fai"
	print >> sys.stderr, "[%s] Convert SAM reference to BAM", timestamp()
	#subprocess.check_output(["samtools", "view", "-o", tmp_dir + "/headered.bam", "-bt", fastidx,  tmp_dir + "/ref.sam"])
	
	p = subprocess.Popen(["samtools", "view", "-o", tmp_dir + "/headered.bam", "-bt", fastidx,  tmp_dir + "/ref.sam"], stderr=subprocess.PIPE)
	out, err = p.communicate()
	
	pattern = re.compile('recognized')
	errlines = err.split("\n")

	for line in errlines:
		m = pattern.search(line)	
		if m:
			print "ERROR: GTF does not match reference fasta"
			print "\tAs reported by samtools:"
			print "\t", line
			quit()
	
	subprocess.call(["samtools","sort",tmp_dir + "/headered.bam",tmp_dir + "/ref"])
	subprocess.call(["samtools", "index", tmp_dir + "/ref.bam"])
	subprocess.call(["rm", tmp_dir + "/headered.bam"])
	subprocess.call(["rm", tmp_dir + "/ref.sam"])
	#subprocess.call(["samtools", "view", "-o", "headered.bam", "-bt", fastidx,  tmp_dir + "/ref.sam"])
	#subprocess.call(["samtools","sort","headered.bam","ref"])
	#subprocess.call(["samtools", "index", "ref.bam"])
	#subprocess.call(["rm", "headered.bam"])
	#subprocess.call(["rm", "ref.sam"])
	return(txdict)

def sam_to_bam(samfile_prefix,fastafile,output_dir):
	'''
	From a sam file and fasta, creates a BAM
	'''
	print "[***************] Converting to BAM format"
	subprocess.call(["samtools", "faidx", fastafile])
	fastidx = fastafile + ".fai"

	mycommands = ["samtools", "view", "-o", output_dir + "/headered.bam", "-bt", fastidx,  output_dir + "/" + samfile_prefix + ".sam"]
	print "[running]", " ".join(mycommands)
	subprocess.call(mycommands)

	mycommands = ["samtools","sort",output_dir + "/headered.bam",output_dir + "/" + samfile_prefix]
	print "[running]", " ".join(mycommands)
	subprocess.call(mycommands)
	
	mycommands = ["samtools", "index", output_dir + "/" + samfile_prefix + ".bam"]
	print "[running]", " ".join(mycommands)
	subprocess.call(mycommands)
	
	mycommands = ["rm", output_dir + "/headered.bam"]
	print "[running]", " ".join(mycommands)
	subprocess.call(mycommands)

def prepare_output_dir(output_dir):
    logging_dir = output_dir + "/logs/"
    tmp_dir = output_dir + "/tmp/"
    print >> sys.stderr, "[%s] Preparing output location %s" % (timestamp(), output_dir)
    if os.path.exists(output_dir):
        pass
    else:        
        os.mkdir(output_dir)
        
    if os.path.exists(logging_dir):
        pass
    else:        
        os.mkdir(logging_dir)
        
    if os.path.exists(tmp_dir):
        pass
    else:        
        os.mkdir(tmp_dir)

def prepare_basic_output_dir(output_dir):
	'''
	Doesn't make a tmp or log directory
	'''
	print >> sys.stderr, "[%s] Preparing output location %s" % (timestamp(), output_dir)
	if os.path.exists(output_dir):
		pass
	else:        
		os.mkdir(output_dir)

if __name__ == "__main__":
    sys.exit(main())







