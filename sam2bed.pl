#!/usr/bin/perl -w

use strict;
use warnings;

use Data::Dumper;
use Carp;
use File::Basename;
use Getopt::Long;


my $prog = basename ($0);
my $pairedEnd = 0;
my $printUniqOnly = 0;
my $verbose = 0;
my $useRNAStrand = 0; # use the strand of the RNA instead of the read

GetOptions (
	"uniq"=>\$printUniqOnly,
	"use-RNA-strand"=>\$useRNAStrand,    
	"v|verbose"=>\$verbose);

if (@ARGV != 2 && @ARGV != 3)
{
	print "Converts OLego SAM format to BED format, works for paired end data and saves into a single BED file, only reports the major alignments. \n";
	print "Usage: $prog [options] <in.sam> <out.bed>\n";
	#print " -p: paired-end data\n";
	print "--uniq : print uniquely mapped reads only\n";
	print "--use-RNA-strand: force to use the strand of the RNA based on the XS tag \n";
	print " -v    : verbose\n";
	exit (1);
}

my ($inSAMFile, $outBedFile) = @ARGV;
#my $outBedFile2 = "";
#if (@ARGV == 3)
#{
#	$outBedFile2 = $ARGV[2];
#	$pairedEnd = 1;
#}

my ($fin, $fout);

open ($fin, "<$inSAMFile") || Carp::croak "cannot open file $inSAMFile to read\n";
open ($fout, ">$outBedFile") || Carp::croak "cannot open file $outBedFile to write\n";
#if ($pairedEnd)
#{
#	open ($fout2, ">$outBedFile2") || Carp::croak "cannot open file $outBedFile2 to write\n";
#}


my $i = 0;
my $found = 0;

while (my $line = <$fin>)
{
	chomp $line;

	next if $line=~/^\s*$/;
	next if $line=~/^\@/;

	print "$i ...\n" if $verbose && $i % 5000 == 0;
	$i++;

	my ($QNAME, $FLAG, $RNAME, $POS, $MAPQ, $CIGAR, $MRNM, $MPOS, $ISIZE, $SEQ, $QUAL, $TAG) = split (/\s+/, $line, 12);

	next if $CIGAR eq '*'; #no alignment

	#print $line, "\n";
	my $flagInfo = decodeSAMFlag ($FLAG);
	#Carp::croak Dumper ($flagInfo), "\n";
	#Carp::croak "inconsistency in specifying PE or SE data\n" if ($flagInfo->{'PE'} + $pairedEnd == 1);
	#next unless $flagInfo->{'query_map'};

	my $strand = $flagInfo->{'query_strand'};
	if($useRNAStrand)
	{
	    if ($TAG=~/XS\:\S*\:([-+\.])/)
	    {
		$strand = $1;
		$strand = '+' if ($1 eq '.');
	    }
	}
	my $read1_or_2 = $flagInfo->{'read_1_or_2'};

	my $name = $QNAME; #substr ($QNAME, 1);
	my $chrom = $RNAME;
	my $chromStart = $POS - 1;

	$TAG = "" unless $TAG;
	my $score = 0;
	if ($TAG=~/NM\:\S*\:(\d+)/)
	{
			#Carp::croak "OK\n";
		$score = $1;
	}

	my $uniq = 0;
	$uniq = 1 if $TAG=~/XT:A:U/;
	

	my ($chromEnd, $block1End, $block2Start, $blockNum);
	my $outStr;

=obsolete
	$TAG=~/\:(\d+)$/;
	my $score = $1;
	
	$blockNum = 1;
	if ($CIGAR =~/^(\d+)M$/)
	{
		$chromEnd = $chromStart + $1 - 1;
		$outStr = join ("\t", $chrom, $chromStart, $chromEnd + 1, $name, $score, $strand);
	}
	elsif ($CIGAR =~/^(\d+)M(\d+)N(\d+)M$/) #two block
	{
		#Carp::croak "junctions detected: $line\n";
		$blockNum = 2;
		$block1End = $chromStart + $1 - 1;
		$block2Start = $chromStart + $1 + $2;
		$chromEnd = $block2Start + $3 - 1;
		$outStr = join ("\t", $chrom, $chromStart, $chromEnd + 1, $name, $score, $strand, 
					$chromStart, $chromEnd + 1, 0, 2, "$1,$3", "0,".($1+$2));
	}
	else
	{
		Carp::croak "unexpected CIGAR string: $CIGAR in $QNAME\n";
	}
=cut

	if ($CIGAR=~/[^\d+|M|N|I|D]/g)
	{
		Carp::croak "unexpected CIGAR string: $CIGAR in $QNAME: $SEQ\n";
	}

	my (@blockSizes, @blockStarts);

	my $currLen = 0;
	my $extendBlock = 0;
	while ($CIGAR=~/(\d+)([M|N|I|D])/g)
	{
		my ($size, $type) = ($1, $2);
		
		if ($type eq 'I' || $type eq 'D')
		{	#insertion in reads
			$extendBlock = 1;
			if ($type eq 'D')
			{
				my $n = @blockSizes;
				if ($n < 1)
				{
					$chromStart += $size;
				}
				else
				{
					#Carp::croak $line, "\n" if @blockSizes <= 0;
					$blockSizes[$#blockSizes] += $size; # if $type eq 'D';
					$currLen += $size;
				}
			}
			next;
		}
		if ($type eq 'M')
		{
			if ($extendBlock && @blockSizes > 0)
			{
				#extend the previous block
				my $n = @blockSizes;
				#Carp::croak $line, "\n" if $n <= 0;
				$blockSizes[$n-1] += $size;
			}
			else
			{
				push @blockSizes, $size;
				push @blockStarts, $currLen;
			}
			$extendBlock = 0;
		}
		$currLen += $size;
	}

	my $blockCount = @blockSizes;
	$chromEnd = $chromStart + $blockStarts[$blockCount-1] + $blockSizes[$blockCount-1] - 1;

	$outStr = join ("\t", $chrom, $chromStart, $chromEnd + 1, $name, $score, $strand,
		$chromStart, $chromEnd + 1, 0, $blockCount, join (",", @blockSizes), join(",", @blockStarts));

#	if ($pairedEnd && $read1_or_2 == 2)
#	{
#		print $fout2 $outStr, "\n" unless $printUniqOnly == 1 && $uniq == 0; # unless $flagInfo->{'mate_map'};
#	}
#	else
#	{
		if ($printUniqOnly == 0 || $uniq == 1)
		{
			print $fout $outStr, "\n" unless $flagInfo->{'query_map'};
		}
#	}
}


close ($fin);
close ($fout);
#close ($fout2) if $pairedEnd;


sub decodeSAMFlag
{
	my $flag = $_[0];
	
	#print "flag = $flag\n";
	$flag = sprintf ("%012b", $flag);

	#print "flag binary = $flag\n";
	my @flags = split (//, $flag);

	my $flagInfo = {
		PE=>$flags[11],
		PE_map=>$flags[10],
		query_map=>$flags[9],
		mate_map=>$flags[8],
		query_strand=>$flags[7] == 0 ? '+' : '-',
		mate_strand=>$flags[6] == 0 ? '+' : '-',
		read_1_or_2=> $flags[5] == 1 ? 1 : 2 };
	return $flagInfo;
}


