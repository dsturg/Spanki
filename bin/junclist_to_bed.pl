#!/usr/bin/perl
#

use strict;
use FileHandle;
use Getopt::Long;

# This version just makes generic anchor 

my $FILENAME= $ARGV[0];
my $anchor = $ARGV[1];

if ($anchor < 1) {
	$anchor = 10
}

#my $READLEN = $ARGV[1]; 
#SEGMENT_ALIGNMENT output file
#my $BED = new FileHandle ">$FILENAME".".bed" or die"can't open $FILENAME".".bed"; #OUTPUT: BROWSER TRAK FILE 
#my $BED_UCSC = new FileHandle ">$FILENAME".".ucsc.bed" or die"can't open $FILENAME".".bed"; #OUTPUT: BROWSER TRAK FILE 
my @f;
my $line;
my $chromo;
my $start;
my $end;
my $strand;

my $juncid;
my $coverage;

my $offset1;
my $offset2;

my $blocksize1;
my $blocksize2;

my $blockstart1;
my $blockstart2;

my @g;
my @h;
#my $readkey;

#my %readhash;

my $leftend;
my $rightstart;

my $intronstart;
my $intronend;
#$BED->printf("track name='$FILENAME'\n"); 
#$BED_UCSC->printf("track name='$FILENAME'\n");     
open(SEQ, "$FILENAME") or die("can't open $FILENAME"); 


while ($line = <SEQ>){    
  
@f = split /\t+/,$line;

	if ($line !~ /juncid/) {
		chomp $line;
		@f = split /\t+/,$line;
		@g = split /[\:]+/,$f[0];
		@h = split /[\_]+/,$g[1];
		$chromo = $g[0];
		#$chromo =~ s/chr//;
		$intronstart = $h[0] - 1;
		$intronend = $h[1];
		$strand = $g[2];
		
		$juncid = $f[0];
		$coverage = $f[1];
		#if ($f[1]) {
		#	$coverage = $f[1];
		#} else {
	#		$coverage = 1
		#}
		
		#print $line,"\n";
		#print "$chromo\t$start\t$leftend\t$rightstart\t$end\n";
		#sleep 2;
		
		$chromo =~ s/>//;
			
		$blocksize1 = $anchor;
		$blocksize2 = $anchor;
		$offset1 = $rightstart - $start;
	    #chrXHet	800	1767	JUNC00000001	2	+	800	1767	255,0,0	2	20,63	0,904
	
		my $intronsize = $intronend - $intronstart;
	
		print "$chromo\t", $intronstart - $anchor, "\t", $intronend + $anchor, "\t", $juncid, "\t", $coverage, "\t";
		print "$strand\t", $intronstart - $anchor, "\t", $intronend + $anchor, "\t255,0,0\t2\t",$anchor,",",$anchor,"\t";
		print "0,", $intronsize + $anchor, "\n";
	
	
	}
  
}
  




