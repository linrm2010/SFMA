#! /usr/bin/perl -w
################################################################################
# split_circular_genome_sequence.pl version 1.0
# Copyright (C) Runmao Lin, Tong Liu, Xiaoting Wang, Fanxing Yang, Zhiyin Wang, 2023
# Contact (E-mail): linrunmao@hainanu.cn
# 
# This program is provided under the terms of a personal license to the recipient and 
# may only be used for the recipient's own research at an academic insititution.
# 
# For using this program in a company or for commercial purposes, a commercial license 
# is required.
# 
# The purpose of this program is to split the circular genome sequence.
################################################################################

use strict;
use warnings;
use Getopt::Long;
use FindBin qw($Bin);
use Data::Dumper;

my $genome;
my $position_site;
my $output;
my $help;
my $version=sprintf("%.1f",1.0);

############################## genome format ###################################
#>mitogenome
#GGCCAAGTTTGATGCATATAATTACCGCAAGTTGCTCCTTGAACATATAA
#ATCGGATCCAGTCATGGTGAAAGCCGCCTTAAGGCCCTCGCCTGGCCGAT
################################################################################

########################### position_of_split_site.txt #########################
########## the 100th base will be the first base in the output file ############
#mitogenomeON764439.1	100
################################################################################

GetOptions
(
	"genome=s" => \$genome,                      # string
	"position_site=s" => \$position_site,        # string
	"output=s" => \$output,                      # string
	"help" => \$help                             # flag
);

################################################################################
# usage
############################# usage begin ######################################
my $usage= << "USAGE";

Example: $0 -genome  mitogenome.ori.fa  -position_site  position_of_split_site.txt  -output  mitogenome.update.fa 
version: 1.0
Options:
	-genome <file>             one file contains genome sequences, such as 'mitogenome.ori.fa'
	-position_site <file>      one file contains position for spliting genome sequences, such as 'position_of_split_site.txt'
	-output <file>             output file, such as 'mitogenome.update.fa'
	-help                      print help information

USAGE
############################## usage end #######################################

if ($help || !(defined $genome) || !(defined $position_site) || !(defined $output))
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

my %seq;
my %title;
my $seq_id="";
open S,"$genome" || die "Cannot open the file '$genome'.\n";
while(<S>)
{
	chomp;
	if($_=~/>/)
	{
		$seq_id=(split(/\s+/,$_))[0];
		$seq_id=~s/>//;
		$title{$seq_id}=$_;
	}
	elsif($_ ne "")
	{
		my $seq_line=$_;
		$seq_line=~s/\s+//g;
		$seq_line=uc $seq_line;
		$seq{$seq_id}.=$seq_line;
	}
}
close S;

print "Reading '$position_site' file ...\n";

my %split_site;
open B,"$position_site" || die "Cannot open the file '$position_site'.\n";
while(<B>)
{
	chomp;
	my @sp=split(/\t/,$_);
	if($_ ne "" && !(defined $split_site{$sp[0]}))
	{
		$split_site{$sp[0]}=$sp[1];
	}
	elsif($_ ne "")
	{
		print "Error! Repetitive occurrence of split site for '$sp[0]'.\n";
		exit;
	}
}
close B;

open OUT,">$output";
foreach my $id(sort keys %seq)
{
	if(defined $split_site{$id})
	{
		my $seq_length=length $seq{$id};
		my $change_seq=substr($seq{$id},$split_site{$id}-1,$seq_length-$split_site{$id}+1);
		$change_seq.=substr($seq{$id},0,$split_site{$id}-1);
		$seq{$id}=$change_seq;
	}
	my $output_seq="";
	my @sp=split(//,$seq{$id});
	for(my $i=0;$i<@sp;$i++)
	{
		$output_seq.=$sp[$i];
		$output_seq.="\n" if(($i+1)%100==0);
	}
	$output_seq=~s/\n$//;
	print OUT "$title{$id}\n";
	print OUT "$output_seq\n";
}
close OUT;

($timesecond, $timeminute, $timehour, $timedaymonths, $timemonth, $timeyear, $timedayweek, $timedayYear, $timedaylightsavings) = localtime();
$timeyear+=1900;
print "[$timehour\:$timeminute\:$timesecond\, $months[$timemonth] $timedaymonths\, $timeyear] ";
print "End of running the \'$0\' program.\n";

__END__
