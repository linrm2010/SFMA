#! /usr/bin/perl -w
################################################################################
# change_degenerate_bases.pl version 1.0
# Copyright (C) Runmao Lin, Tong Liu, Xiaoting Wang, Fanxing Yang, Zhiyin Wang, 2023
# Contact (E-mail): linrunmao@hainanu.cn
# 
# This program is provided under the terms of a personal license to the recipient and 
# may only be used for the recipient's own research at an academic insititution.
# 
# For using this program in a company or for commercial purposes, a commercial license 
# is required.
# 
# The purpose of this program is to identify degenerate bases in the genome sequences.
################################################################################

use strict;
use warnings;
use Getopt::Long;
use FindBin qw($Bin);
use Data::Dumper;

my $genome;
my $output;
my $help;
my $version=sprintf("%.1f",1.0);

############################## genome format ###################################
#>mitogenome
#GGCCAAGTTTGATGCATATAATTACCGCAAGTTGCTCCTTGAACATATAA
#ATCGGATCCAGTCATGGTGAAAGCCGCCTTAAGGCCCTCGCCTGGCCGAT
################################################################################

GetOptions
(
	"genome=s" => \$genome,              # string
	"output=s" => \$output,              # string
	"help" => \$help                     # flag
);

################################################################################
# usage
############################# usage begin ######################################
my $usage= << "USAGE";

Example: $0 -genome  genome.ori.fa  -output  genome.update.fa 
version: 1.0
Options:
	-genome <file>      one file contains genome sequences, such as 'mitogenome.ori.fa'
	-output <file>      output file, such as 'genome.update.fa'
	-help               print help information

USAGE
############################## usage end #######################################

if ($help || !(defined $genome) || !(defined $output))
{
	print $usage;
	exit;
}

my @months = qw(Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec);
my ($timesecond, $timeminute, $timehour, $timedaymonths, $timemonth, $timeyear, $timedayweek, $timedayYear, $timedaylightsavings) = localtime();
$timeyear+=1900;

print "[$timehour\:$timeminute\:$timesecond\, $months[$timemonth] $timedaymonths\, $timeyear] ";
print "Start to run the \'$0\' program ...\n";

################################################################################
# Reading 'genome'
################################################################################

################################ degenerate_bases ##############################
# A+T	W
# C+G	S
# T+G	K
# A+C	M
# C+T	Y
# A+G	R
# A+C+G	V
# A+T+G	D
# T+C+G	B
# A+T+C	H
# A+G+C+T	N
################################################################################

print "Reading '$genome' file ...\n";

open OUT,">$output";

my $seq_id="";
my $degenerate_bases="";
my %length;
open S,"$genome" || die "Cannot open the file '$genome'.\n";
while(<S>)
{
	chomp;
	if($_=~/>/)
	{
		print OUT "$_\n";
		$seq_id=(split(/\s+/,$_))[0];
		$seq_id=~s/>//;
	}
	elsif($_ ne "")
	{
		my $seq_line=$_;
		$seq_line=uc $seq_line;
		$seq_line=~s/\s+//g;
		my @sp=split(//,$seq_line);
		for(my $i=0;$i<@sp;$i++)
		{
			$length{$seq_id}++;
			if($sp[$i] eq "W")
			{
				$degenerate_bases.="$seq_id\t$length{$seq_id}\tW\tA\n";
				$sp[$i]="A";
			}
			elsif($sp[$i] eq "S")
			{
				$degenerate_bases.="$seq_id\t$length{$seq_id}\tS\tC\n";
				$sp[$i]="C";
			}
			elsif($sp[$i] eq "K")
			{
				$degenerate_bases.="$seq_id\t$length{$seq_id}\tK\tT\n";
				$sp[$i]="T";
			}
			elsif($sp[$i] eq "M")
			{
				$degenerate_bases.="$seq_id\t$length{$seq_id}\tM\tA\n";
				$sp[$i]="A";
			}
			elsif($sp[$i] eq "Y")
			{
				$degenerate_bases.="$seq_id\t$length{$seq_id}\tY\tC\n";
				$sp[$i]="C";
			}
			elsif($sp[$i] eq "R")
			{
				$degenerate_bases.="$seq_id\t$length{$seq_id}\tR\tA\n";
				$sp[$i]="A";
			}
			elsif($sp[$i] eq "V")
			{
				$degenerate_bases.="$seq_id\t$length{$seq_id}\tV\tA\n";
				$sp[$i]="A";
			}
			elsif($sp[$i] eq "D")
			{
				$degenerate_bases.="$seq_id\t$length{$seq_id}\tD\tA\n";
				$sp[$i]="A";
			}
			elsif($sp[$i] eq "B")
			{
				$degenerate_bases.="$seq_id\t$length{$seq_id}\tB\tT\n";
				$sp[$i]="T";
			}
			elsif($sp[$i] eq "H")
			{
				$degenerate_bases.="$seq_id\t$length{$seq_id}\tH\tA\n";
				$sp[$i]="A";
			}
		}
		$seq_line=join("",@sp);
		print OUT "$seq_line\n";
	}
}
close S;
close OUT;

if($degenerate_bases ne "")
{
	open OB,">$output\.degenerate_bases.txt";
	print OB "SequenceID\tPosition\tDegenerateBase\tChangedBase\n";
	print OB $degenerate_bases;
	close OB;
}

($timesecond, $timeminute, $timehour, $timedaymonths, $timemonth, $timeyear, $timedayweek, $timedayYear, $timedaylightsavings) = localtime();
$timeyear+=1900;
print "[$timehour\:$timeminute\:$timesecond\, $months[$timemonth] $timedaymonths\, $timeyear] ";
print "End of running the \'$0\' program.\n";

__END__
