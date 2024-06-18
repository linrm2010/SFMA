#! /usr/bin/perl -w
################################################################################
# reverse_complementary.pl version 1.0
# Copyright (C) Runmao Lin, Tong Liu, Xiaoting Wang, Fanxing Yang, Zhiyin Wang, 2023
# Contact (E-mail): linrunmao@hainanu.cn
# 
# This program is provided under the terms of a personal license to the recipient and 
# may only be used for the recipient's own research at an academic insititution.
# 
# For using this program in a company or for commercial purposes, a commercial license 
# is required.
# 
# The purpose of this program is to obtain the reverse and complementary sequences of genome.
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
	"genome=s" => \$genome,                      # string
	"output=s" => \$output,                      # string
	"help" => \$help                             # flag
);

################################################################################
# usage
############################# usage begin ######################################
my $usage= << "USAGE";

Example: $0 -genome  mitogenome.ori.fa  -output  mitogenome.rev_com.fa 
version: 1.0
Options:
	-genome <file>             one file contains genome sequences, such as 'mitogenome.ori.fa'
	-output <file>             the reverse and complementary sequences of genome, such as 'mitogenome.rev_com.fa'
	-help                      print help information

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

print "Reading '$genome' file ...\n";

my $title='';
my ($seq,$rev_com_seq)=();
open IN,"$genome" || die "Cannot open the file '$genome'.\n";
while(<IN>)
{
	chomp;
	if($_=~/^\>/)
	{
		$title=$_;
	}
	else
	{
		$seq->{$title}.=$_;
	}
}
close IN;
$rev_com_seq=$seq;
open OUT,">$output";
foreach $title(keys %$rev_com_seq)
{
	$rev_com_seq->{$title}=uc($rev_com_seq->{$title});
	$rev_com_seq->{$title}=~s/A/K/g;
	$rev_com_seq->{$title}=~s/C/L/g;
	$rev_com_seq->{$title}=~s/T/A/g;
	$rev_com_seq->{$title}=~s/G/C/g;
	$rev_com_seq->{$title}=~s/K/T/g;
	$rev_com_seq->{$title}=~s/L/G/g;
	my @temp_seq=split(//,$rev_com_seq->{$title});
	@temp_seq=reverse(@temp_seq);
#	$rev_com_seq->{$title}=join('',@temp_seq);
	my $rev_seq='';
	for(my $i=0;$i<@temp_seq;$i++)
	{
		$rev_seq.=$temp_seq[$i];
		if(($i+1)%50==0)
		{
			$rev_seq.="\n";
		}
	}
	print OUT "$title\n$rev_seq";
}
close OUT;

($timesecond, $timeminute, $timehour, $timedaymonths, $timemonth, $timeyear, $timedayweek, $timedayYear, $timedaylightsavings) = localtime();
$timeyear+=1900;
print "[$timehour\:$timeminute\:$timesecond\, $months[$timemonth] $timedaymonths\, $timeyear] ";
print "End of running the \'$0\' program.\n";

__END__
