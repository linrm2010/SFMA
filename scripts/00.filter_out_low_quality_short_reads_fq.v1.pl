#! /usr/bin/perl -w
################################################################################
# filter_out_low_quality_short_reads_fq.pl version 1.0
# Copyright (C) Runmao Lin, Tong Liu, Xiaoting Wang, Fanxing Yang, Zhiyin Wang, 2023
# Contact (E-mail): linrunmao@hainanu.cn
# 
# This program is provided under the terms of a personal license to the recipient and 
# may only be used for the recipient's own research at an academic insititution.
# 
# For using this program in a company or for commercial purposes, a commercial license 
# is required.
# 
# The purpose of this program is to filter out short reads with low quality.
################################################################################

use strict;
use warnings;
use Getopt::Long;
use FindBin qw($Bin);
use Data::Dumper;

my $fq_file;
my $output;
my $min_base_length;
my $min_quality;
my $phred_score;
my $help;
my $version=sprintf("%.1f",1.0);

############################## fq_file format ##################################
#@A00881:406:HHGTHDSXY:3:1101:3269:1297 1:N:0:ATCCTTGG+TAGCCACT
#GGCCAAGTTTGATGCATATAATTACCGCAAGTTGCTCCTTGAACATATAAATCGGATCCAGTCATGGTGAAAGCCGCCTTAAGGCCCTCGCCTGGCCGAT
#+
#FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
################################################################################

GetOptions
(
	"fq_file=s" => \$fq_file,                          # string
	"min_base_length=f" => \$min_base_length,          # numeric
	"min_quality=f" => \$min_quality,                  # numeric
	"phred_score=f"=> \$phred_score,                   # numeric
	"output=s" => \$output,                            # string
	"help" => \$help                                   # flag
);

################################################################################
# usage
############################# usage begin ######################################
my $usage= << "USAGE";

Example: $0 -fq_file  sample_1.fq.gz  -min_base_length  100  -min_quality  20  -phred_score  33  -output  sample_1.filter.fq 
version: 1.0
Options:
	-fq_file <file>               one file contains sequenced data, such as 'sample_1.fq.gz'
	-min_base_length <numeric>    such as '100'
	-min_quality <numeric>        such as '20'
	-phred_score <numeric>        such as '33'
	-output <file>                output file, such as 'sample_1.filter.fq'
	-help                         print help information

USAGE
############################## usage end #######################################

if ($help || !(defined $fq_file) || !(defined $min_base_length) || !(defined $min_quality) || !(defined $phred_score) || !(defined $output))
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
# Reading 'fq_file'
################################################################################

print "Reading the '$fq_file' file ...\n";

open OUT,">$output";
if($fq_file=~/\.gz$/)
{
	open FQ,"gzip -dc $fq_file | " || die "Cannot open the file '$fq_file'.\n";
}
else
{
	open FQ,"$fq_file" || die "Cannot open the file '$fq_file'.\n";
}
while(<FQ>)
{
	chomp;
	if($_=~/^\@/)
	{
		my $name=$_;
		my $seq=<FQ>;
		chomp $seq;
		my @spseq=split(//,$seq);
		my $tag=0;
		my $bpos=0;
		for(my $i=0;$i<@spseq;$i++)
		{
			if($spseq[$i] ne "N" && $tag==0)
			{
				$bpos=$i;
				$tag++;
			}
			elsif($spseq[$i] ne "N" && $tag>0)
			{
				$tag++;
			}
			elsif($spseq[$i] eq "N" && $tag<$min_base_length)
			{
				$tag=0;
				$bpos=0;
			}
			elsif($spseq[$i] eq "N" && $tag>=$min_base_length)
			{
				last;
			}
		}
		if($tag>=$min_base_length)
		{
			$seq=substr($seq,$bpos,$tag);
			<FQ>;
			my $qual=<FQ>;
			chomp $qual;
			$qual=substr($qual,$bpos,$tag);
			my @spqual=split(//,$qual);
			my $qtag=0;
			my $qbpos=0;
			for(my $s=0;$s<@spqual;$s++)
			{
				my $value=ord($spqual[$s])-$phred_score;
				if($value>=20 && $qtag==0)
				{
					$qbpos=$s;
					$qtag++;
				}
				elsif($value>=20 && $qtag>0)
				{
					$qtag++;
				}
				elsif($value<20 && $qtag<$min_base_length)
				{
					$qbpos=0;
					$qtag=0;
				}
				elsif($value<20 && $qtag>=$min_base_length)
				{
					last;
				}
			}
			if($qtag>=$min_base_length)
			{
				$seq=substr($seq,$qbpos,$qtag);
				$qual=substr($qual,$qbpos,$qtag);
				print OUT "$name\n$seq\n+\n$qual\n";
			}
		}
		else
		{
			<FQ>;
			<FQ>;
		}
	}
}
close FQ;
close OUT;

($timesecond, $timeminute, $timehour, $timedaymonths, $timemonth, $timeyear, $timedayweek, $timedayYear, $timedaylightsavings) = localtime();
$timeyear+=1900;
print "[$timehour\:$timeminute\:$timesecond\, $months[$timemonth] $timedaymonths\, $timeyear] ";
print "End of running the \'$0\' program.\n";

__END__
