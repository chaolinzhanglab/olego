
use strict;
use warnings;

use File::Basename;
use Getopt::Long;

my $prog = basename ($0);

my $verbose = 0;
my $strandaware = 0;

GetOptions (
    "v|verbose"=>\$verbose,
    "s"=>\$strandaware
);

if (@ARGV != 2)
{
    print "bed format to junctions\n";
    print "usage: $prog [options] <in.bed> <out> \n";
    exit(1);
}

my ($inputfilename, $outputfilename) = @ARGV;

my ($fin, $fout);

my %junchash;

open($fin, $inputfilename) or Carp::croak "cannot open file $inputfilename to read!\n";
while(my $line = <$fin>)
{
    chomp($line);
    my @a =split("\t", $line);
    my @blockSizes = split(",", $a[10]);
    my @blockStarts = split(",", $a[11]);

    for (my $i=0; $i<$a[9]-1; $i++)
    {
	    my $start = $a[1] + $blockStarts[$i] + $blockSizes[$i];
	    my $end = ($a[1] + $blockStarts[$i+1]);
	    if( not $strandaware){
		$junchash{join(",", $a[0], $start, $end)}++;
	    }
	    else
	    {
		$junchash{join(",", $a[0], $start, $end, $a[5])}++;
	    }
	    
    }

}
close($fin);

open($fout, ">$outputfilename") or Carp::croak "cannot open file $outputfilename to write!\n";
foreach my $key (keys %junchash)
{
    my @a = split(",", $key);
    print $fout $a[0],"\t";
    print $fout $a[1], "\t";
    print $fout $a[2], "\t";
    print $fout "\.\t";
    print $fout $junchash {$key},"\t";
    if ( not $strandaware){
	print $fout "+\n";
    }
    else
    {
	print $fout $a[3],"\n";
    }
    #print $a[5],"\n";
}
close($fout);
