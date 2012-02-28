#!/usr/bin/perl -w

use strict;
use warnings;

use Getopt::Long;
use Carp;


my $verbose = 0;

my $maxdist = 5000000; #cross 11 exons

GetOptions (
    "d|distance"=>\$maxdist,
    "v|verbose"=>\$verbose
);

my $prog = $0;

if ( @ARGV != 3)
{
    print "Merge sam format output from PE reads\n";
    print "Usage: $prog [options] <end1.sam> <end2.sam> <out.sam>\n";
    print "-d	max distance between the two ends on genome [$maxdist]\n";
    print "-v	verbose\n";
    exit 1;
}


my ($inSAMFile1,$inSAMFile2, $outFile) = @ARGV;
my ($fin1,$fin2, $fout); 

open ($fin1, "<$inSAMFile1") || Carp::croak "cannot open file $inSAMFile1 to read\n";
open ($fin2, "<$inSAMFile2") || Carp::croak "cannot open file $inSAMFile2 to read\n";

open ($fout, ">$outFile") || Carp::croak "cannot open file $outFile to write\n";


while (my $line1 = <$fin1>)
{
    chomp $line1;
    my $line2 = <$fin2>;
    chomp $line2;
    if ($line1=~/^\@/)
    {
	print $fout $line1,"\n";
	next;
    }
    #get information from file 1
    my ($QNAME1, $FLAG1, $RNAME1, $POS1, $MAPQ1, $CIGAR1, $MRNM1, $MPOS1, $ISIZE1, $SEQ1, $QUAL1, $TAG1) = split (/\s+/, $line1, 12);

    my ($QNAME2, $FLAG2, $RNAME2, $POS2, $MAPQ2, $CIGAR2, $MRNM2, $MPOS2, $ISIZE2, $SEQ2, $QUAL2, $TAG2) = split (/\s+/, $line2, 12);

    # put the infor into the arrays 
    my (@chr1, @chr2);
    my (@pos1, @pos2);
    my (@tag1, @tag2);
    my (@cigar1, @cigar2);
    my (@nm1, @nm2);
    my (@strand1, @strand2);

    push (@chr1, $RNAME1);
    push (@chr2, $RNAME2);
    
    push (@pos1, $POS1);
    push (@pos2, $POS2);

    push (@cigar1, $CIGAR1);
    push (@cigar2, $CIGAR2);
#   modify FLAGs

    $FLAG1 = $FLAG1 | 0x0001 ; #the read is paired
    $FLAG2 = $FLAG2 | 0x0001 ;
    $FLAG1 = $FLAG1 | 0x0040 ; # the read is the first read in the pair
    $FLAG2 = $FLAG2 | 0x0080 ; # the read is the second read in the pair
    
    #work on the flags later after changing the main code
    
    my $flaginfo1 = decodeSAMFlag ($FLAG1);
    my $flaginfo2 = decodeSAMFlag ($FLAG2);
 
    # do not look at strand now
    unless ($TAG1 and $TAG2)
    {
	print $fout $line1,"\n";
	print $fout $line2, "\n";
	next;
    }
    
    push (@strand1, $flaginfo1->{'query_strand'});
     push (@strand2, $flaginfo2->{'query_strand'});
     
    push (@tag1, $TAG1); 

    if( $TAG1=~/NM\:\S*\:(\d+)/)
    {
	push(@nm1, $1);
    }
    else
    {
	push(@nm1, -1);
    }
    
    push (@tag2, $TAG2);

    if( $TAG2=~/NM\:\S*\:(\d+)/)
    {
	push(@nm2, $1);
    }
    else
    {
	push(@nm2, -1);	
    }
    
    # scan XA tags
    # XA should be the last tag

    if ($TAG1=~/XA\:\S\:(\S*)$/)
    {
	my @XAstrs = split(";",$1);
	for(my $i=0; $i<@XAstrs; $i++)
	{
#	    XA:Z:chr4,+149621574,100M,0;chr2,+80678177,100M,0;
	    $XAstrs[$i]=~/^(\w*?),([\+\-])(\d*?),(\w*?),(\d*?)$/;
	    push(@chr1, $1);
	    push(@strand1, $2);
	    push(@pos1, $3);
	    push(@cigar1, $4);
	    push(@nm1, $5);
	}
    }
    if ($TAG2=~/XA\:\S\:(\S*)$/)
      {
	  my @XAstrs = split(";",$1);
	  for(my $i=0; $i<@XAstrs; $i++)
	  {
#	    XA:Z:chr4,+149621574,100M,0;chr2,+80678177,100M,0;
	      $XAstrs[$i]=~/^(\w*?),([\+\-])(\d*?),(\w*?),(\d*?)$/;
	      push(@chr2, $1);
	      push(@strand2, $2);
	      push(@pos2, $3);
	      push(@cigar2, $4);
	      push(@nm2, $5);
	  }
      }
    # scan for the best match
#    my $foundmatches=0, $bestdistance = $distance, $bestnm = 0;
    my %distanceij;
    my %rankij;
    for (my $i = 0; $i<@chr1; $i++)
     {
	 for(my $j = 0; $j<@chr2; $j++)
	{
	    if($chr1[$i] eq $chr2[$j] and abs($pos2[$j]-$pos1[$i])<$maxdist)
	    {
		#$bestdistance =  abs($pos2[$j]-$pos2[$i]);
		#$bestnm = $nm1[$i] + $nm2[$j];
		$distanceij{$i.",".$j} = abs($pos2[$j]-$pos1[$i]);
		$rankij{$i.",".$j} = $i+$j;
		#print join("\t", $i,$j),"\n";	
	    }
	}
     }
     
     if(scalar (keys %distanceij) >0 )
    {
	# only output reliable hits if matches found
	# but some entries could be output several times. 
	my $ontop = 1;
	my $outputline1;
	my $outputline2;
	
	foreach my $ij (sort {$rankij{$a} <=> $rankij{$b}} keys %rankij)
	{
	    my ($i,$j) = split(",", $ij);
	    if($ontop == 1)
	    {
		$outputline1 = join("\t", $QNAME1, $FLAG1, $chr1[$i], $pos1[$i], $MAPQ1, $cigar1[$i], $MRNM1, $MPOS1, $ISIZE1, $SEQ1, $QUAL1, "NM:i:".$nm1[$i]);
		$outputline1 = $outputline1."\tXA:Z:" if(scalar (keys %distanceij) >1);
		$outputline2 = join("\t", $QNAME2, $FLAG2, $chr2[$j], $pos2[$j], $MAPQ2, $cigar2[$j], $MRNM2, $MPOS2, $ISIZE2, $SEQ2, $QUAL2, "NM:i:".$nm2[$j]);
		$outputline2= $outputline2."\tXA:Z:" if(scalar (keys %distanceij) >1);
		
		$ontop = 0;
	    }
	    else
	    {
		$outputline1 = $outputline1 .join(",",$chr1[$i], $strand1[$i].$pos1[$i],$cigar1[$i], $nm1[$i] ).";";
		$outputline2 = $outputline2 .join(",",$chr2[$j], $strand2[$j].$pos2[$j],$cigar2[$j], $nm2[$j] ).";";
	    }
	    
        }
	print $fout $outputline1, "\n";
	print $fout $outputline2, "\n";

    }
    else
    {
	# no match, output the original hits
	print $fout $line1,"\n";
	print $fout $line2,"\n";
	
    }
}

close($fin1);
close($fin2);


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


