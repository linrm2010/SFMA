#! /usr/bin/perl -w
################################################################################
# sequence_statistic.pl version 1.0
# Copyright (C) Runmao Lin, Tong Liu, Xiaoting Wang, Fanxing Yang, Zhiyin Wang, 2023
# Contact (E-mail): linrunmao@hainanu.cn
#
# This program is provided under the terms of a personal license to the recipient and
# may only be used for the recipient's own research at an academic insititution.
#
# For using this program in a company or for commercial purposes, a commercial license
# is required.
#
# The purpose of this program is to sort sequence by length, and report statistics of sequences.
################################################################################

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;

my $sequence;
my $output_file;
my $help;

GetOptions
(
	"genome=s" => \$genome,                                             # string
	"output_file=s" => \$output_file,                                   # string
	"help" => \$help                                                    # flag
);

################################################################################
# usage
############################# usage begin ######################################
my $usage= << "USAGE";

Example: perl $0 -sequence  genome.fa  -output_file  genome.fa.stat 
version: 1.0
Options:
        -sequence <file>                         genome sequences, such as 'genome.fa'
        -output_file <file>                    output file, such as 'genome.fa.stat'
        -help                                  print help information

USAGE
############################## usage end #######################################

if ($help || !(defined $sequence) || !(defined $output_file))
{
	print $usage;
	exit;
}

my @months = qw(Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec);
my ($timesecond, $timeminute, $timehour, $timedaymonths, $timemonth, $timeyear, $timedayweek, $timedayYear, $timedaylightsavings) = localtime();
$timeyear+=1900;

print "[$timehour\:$timeminute\:$timesecond\, $months[$timemonth] $timedaymonths\, $timeyear] ";
print "Start to run the \'$0\' program ...\n";

my %seq;
my %length;
my $slen=0; # sum of length
my $snum=0; # scaffold number
my $snum1k=0;
my $cnum=0; # contig number
my $cnum1k=0;
my $smax=0; # longest scaffold
my $smaxid=""; # id
my $smin=0; # minimal scaffold
my $sminid=""; # id
my $cmax=0; # longest contig
my $cmaxid=""; # id
my $cmin=0; # minimal contig
my $cminid=""; # id
my $n50=0;
my $n90=0;
my $gc=0;
my $gcnum=0;
my $atnum=0;
open IN,"$sequence" || die "Cannot open the file '$sequence'.\n";
$/=">";
<IN>;
$/="\n";
while(<IN>)
{
	my $ti=$_;
	chomp $ti;
	my $id=(split(/\s+/,$ti))[0];
	$/=">";
	$seq{$id}=<IN>;
	chomp $seq{$id};
	$seq{$id}=uc $seq{$id};
	my $tempseq=$seq{$id};
	$tempseq=~s/\n$//;
	$seq{$id}=~s/\s+//g;
	$length{$id}=length $seq{$id};
	my @tcn=split(/N+/,$seq{$id});
	$snum++;
	$snum1k++ if($length{$id}>=1000);
	$slen+=$length{$id};
	if($smax<$length{$id})
	{
		$smax=$length{$id};
		$smaxid=$id;
	}
	if($smin==0 || $smin>$length{$id})
	{
		$smin=$length{$id};
		$sminid=$id;
	}
	for(my $i=0;$i<@tcn;$i++)
	{
		my $tlen=length $tcn[$i];
		$cnum++;
		$cnum1k++ if($tlen>=1000);
		if($cmax<$tlen)
		{
			$cmax=$tlen;
			$cmaxid="$id\-contig".($i+1);
		}
		if($cmin==0 || $cmin>$tlen)
		{
			$cmin=$tlen;
			$cminid="$id\-contig".($i+1);
		}
		my $tgc=$tcn[$i];
		$tgc=~s/[GC]//g;
		my $tatlen=length $tgc;
		$atnum+=$tatlen;
		$gcnum+=$tlen-$tatlen;
	}
	$seq{$id}="\>$ti\n$tempseq";
	$/="\n";
}
close IN;
$gc=$gcnum/($gcnum+$atnum);
my $n50tag=0;
my $n90tag=0;
my $countlen=0;
foreach my $id(sort{$length{$b}<=>$length{$a}} keys %length)
{
	$countlen+=$length{$id};
	if($countlen/$slen>=0.5 && $n50tag==0)
	{
		$n50tag=1;
		$n50=$length{$id};
	}
	if($countlen/$slen>=0.9 && $n90tag==0)
	{
		$n90tag=1;
		$n90=$length{$id};
	}
}
open OB,">$output_file";
printf OB "scaffold num: %10.0f\tscaffold num (>=1000): %8.0f\tscaffold length (bp): %16.0f\n",$snum,$snum1k,$slen;
printf OB "contig num: %12.0f\tcontig num (>=1000): %10.0f\t\n",$cnum,$cnum1k; 
printf OB "max scaffold length (bp): %10.0f\tmax scaffold id: %16s\n",$smax,$smaxid;
printf OB "min scaffold length (bp): %10.0f\tmin scaffold id: %16s\n",$smin,$sminid;
printf OB "scaffold N50 (bp): %10.0f\tscaffold N90 (bp): %10.0f\tGC: %10.4f\n",$n50,$n90,$gc;
close OB;

($timesecond, $timeminute, $timehour, $timedaymonths, $timemonth, $timeyear, $timedayweek, $timedayYear, $timedaylightsavings) = localtime();
$timeyear+=1900;
print "[$timehour\:$timeminute\:$timesecond\, $months[$timemonth] $timedaymonths\, $timeyear] ";
print "End of running the \'$0\' program.\n";

__END__
