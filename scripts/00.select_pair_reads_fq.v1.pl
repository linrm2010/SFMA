#! /usr/bin/perl -w
################################################################################
# select_pair_reads_fq.pl version 1.0
# Copyright (C) Runmao Lin, Tong Liu, Xiaoting Wang, Fanxing Yang, Zhiyin Wang, 2023
# Contact (E-mail): linrunmao@hainanu.cn
# 
# This program is provided under the terms of a personal license to the recipient and 
# may only be used for the recipient's own research at an academic insititution.
# 
# For using this program in a company or for commercial purposes, a commercial license 
# is required.
# 
# The purpose of this program is to identify paired-end reads.
################################################################################

use strict;
use warnings;
use Getopt::Long;
use FindBin qw($Bin);
use Data::Dumper;

my $file_1_fq;
my $file_2_fq;
my $output_prefix;
my $help;
my $version=sprintf("%.1f",1.0);

####################### file_1_fq, file_2_fq format ############################
#@A00881:406:HHGTHDSXY:3:1101:3269:1297 1:N:0:ATCCTTGG+TAGCCACT
#GGCCAAGTTTGATGCATATAATTACCGCAAGTTGCTCCTTGAACATATAAATCGGATCCAGTCATGGTGAAAGCCGCCTTAAGGCCCTCGCCTGGCCGAT
#+
#FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
################################################################################

GetOptions
(
	"file_1_fq=s" => \$file_1_fq,                      # string
	"file_2_fq=s" => \$file_2_fq,                      # string
	"output_prefix=s" => \$output_prefix,              # string
	"help" => \$help                                   # flag
);

################################################################################
# usage
############################# usage begin ######################################
my $usage= << "USAGE";

Example: $0 -file_1_fq  sample_1.filter.fq.gz  -file_2_fq  sample_2.filter.fq.gz  -output_prefix  sample 
version: 1.0
Options:
	-file_1_fq <file>             one file contains sequenced data, such as 'sample_1.filter.fq.gz'
	-file_2_fq <file>             one file contains sequenced data, such as 'sample_2.filter.fq.gz'
	-output_prefix <strings>      output prefix, such as 'sample'
	-help                         print help information

USAGE
############################## usage end #######################################

if ($help || !(defined $file_1_fq) || !(defined $file_2_fq) || !(defined $output_prefix))
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
# Reading 'file_1_fq'
################################################################################

print "Reading '$file_1_fq' file ...\n";

my %reads;
my %pair;
open A,"$file_1_fq" || die "Cannot open the file '$file_1_fq'.\n";
while(<A>)
{
	chomp;
	if($_=~/^\@/)
	{
		my $id=(split(/\s+/,$_))[0];
		if($id=~/\/[12]$/)
		{
			my $seq=<A>;
			chomp $seq;
			<A>;
			my $qual=<A>;
			chomp $qual;
			my $tempid=$id;
			$tempid=~s/\/[12]$//;
			$reads{$tempid}="$id\n$seq\n+\n$qual\n";
			$pair{$tempid}=1;
		}
		else
		{
			my $seq=<A>;
			chomp $seq;
			<A>;
			my $qual=<A>;
			chomp $qual;
			$reads{$id}="$id\/1\n$seq\n+\n$qual\n";
			$pair{$id}=1;
		}
	}
}
close A;

open OA,">$output_prefix\_1.fq";
open OB,">$output_prefix\_2.fq";
open OC,">$output_prefix\_single_1.fq";
open OD,">$output_prefix\_single_2.fq";

################################################################################
# Reading 'file_2_fq'
################################################################################

print "Reading '$file_2_fq' file ...\n";

open B,"$file_2_fq" || die "Cannot open the file '$file_2_fq'.\n";
while(<B>)
{
	chomp;
	if($_=~/^\@/)
	{
		my $id=(split(/\s+/,$_))[0];
		if($id=~/\/[12]$/)
		{
			my $seq=<B>;
			chomp $seq;
			<B>;
			my $qual=<B>;
			chomp $qual;
			my $tempid=$id;
			$tempid=~s/\/[12]$//;
			if(defined $reads{$tempid})
			{
				print OA "$reads{$tempid}";
				print OB "$id\n$seq\n+\n$qual\n";
				$pair{$tempid}++;
			}
			else
			{
				print OD "$id\n$seq\n+\n$qual\n";
			}
		}
		else
		{
			my $seq=<B>;
			chomp $seq;
			<B>;
			my $qual=<B>;
			chomp $qual;
			if(defined $reads{$id})
			{
				print OA "$reads{$id}";
				print OB "$id\/2\n$seq\n+\n$qual\n";
				$pair{$id}++;
			}
			else
			{
				print OD "$id\/2\n$seq\n+\n$qual\n";
			}
		}
	}
}
close B;
foreach my $id(keys %pair)
{
	if($pair{$id}==1)
	{
		print OC "$reads{$id}";
	}
}
close OA;
close OB;
close OC;
close OD;

($timesecond, $timeminute, $timehour, $timedaymonths, $timemonth, $timeyear, $timedayweek, $timedayYear, $timedaylightsavings) = localtime();
$timeyear+=1900;
print "[$timehour\:$timeminute\:$timesecond\, $months[$timemonth] $timedaymonths\, $timeyear] ";
print "End of running the \'$0\' program.\n";

__END__
