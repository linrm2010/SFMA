#! /usr/bin/perl -w
################################################################################
# integrate_palindrome.pl version 1.0
# Copyright (C) Runmao Lin, Tong Liu, Xiaoting Wang, Fanxing Yang, Zhiyin Wang, 2023
# Contact (E-mail): linrunmao@hainanu.cn
#
# This program is provided under the terms of a personal license to the recipient and
# may only be used for the recipient's own research at an academic insititution.
#
# For using this program in a company or for commercial purposes, a commercial license
# is required.
#
# The purpose of this program is to integrate information of palindrome sequences.
################################################################################

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;

my $gene_gff;
my $tRNA_gff;
my $rRNA_gff;
my $DNA_Analyzer_Palindrome;
my $output_file;
my $help;

GetOptions
(
	"gene_gff=s" => \$gene_gff,                                         # string
	"tRNA_gff=s" => \$tRNA_gff,                                         # string
	"rRNA_gff=s" => \$rRNA_gff,                                         # string
	"DNA_Analyzer_Palindrome=s" => \$DNA_Analyzer_Palindrome,           # string
	"output_file=s" => \$output_file,                                   # string
	"help" => \$help                                                    # flag
);

################################################################################
# usage
############################# usage begin ######################################
my $usage= << "USAGE";

Example: perl $0  -gene_gff  integrate.gene.gff  -tRNA_gff  tRNA_result.gff  -rRNA_gff  rRNA.gff  -DNA_Analyzer_Palindrome DNA_Analyzer_Palindrome.txt  -output_file  DNA_Analyzer_Palindrome.check.txt 
version: 1.0
Options:
        -gene_gff <file>                 gene gff file, such as 'integrate.gene.gff'
        -tRNA_gff <file>                 mfanno predicted tRNA, such as 'tRNA_result.gff'
        -rRNA_gff <file>                 tRNAscan predicted tRNA, such as 'rRNA.gff'
        -DNA_Analyzer_Palindrome <file>  tRNAscan predicted tRNA, such as 'DNA_Analyzer_Palindrome.txt'
        -output_file <strings>           output file, such as 'DNA_Analyzer_Palindrome.check.txt'
        -help                            print help information

USAGE
############################## usage end #######################################

if ($help || !(defined $gene_gff) || !(defined $tRNA_gff) || !(defined $rRNA_gff) || !(defined $DNA_Analyzer_Palindrome) || !(defined $output_file))
{
	print $usage;
	exit;
}

my @months = qw(Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec);
my ($timesecond, $timeminute, $timehour, $timedaymonths, $timemonth, $timeyear, $timedayweek, $timedayYear, $timedaylightsavings) = localtime();
$timeyear+=1900;

print "[$timehour\:$timeminute\:$timesecond\, $months[$timemonth] $timedaymonths\, $timeyear] ";
print "Start to run the \'$0\' program ...\n";

my %pos;

print "Reading the '$rRNA_gff' file ...\n";
# CP084950.1	mitfi	rRNA	8149	12824	.	+	.	ID=rRNA:rnl;Name=rnl
# CP084950.1	mitfi	intron	10475	12222	.	+	.	ID=rRNA:rnl;Name=rnl
open RR,"$rRNA_gff" || die "Cannot open the file '$rRNA_gff'.\n";
while(<RR>)
{
	chomp;
	my @sp=split(/\t/,$_);
	if($_=~/\trRNA\t/)
	{
		for(my $i=$sp[3];$i<=$sp[4];$i++)
		{
			$pos{$i}=1;
		}
	}
	if($_=~/\tintron\t/)
	{
		for(my $i=$sp[3];$i<=$sp[4];$i++)
		{
			$pos{$i}=0;
		}
	}
}
close RR;

print "Reading the '$tRNA_gff' file ...\n";
# CP084950.1	tRNAscan-SE	tRNA	1672 	1740 	.	+	.	ID=tRNA:Phe(AAA)_1;Name=Phe(AAA)_1
open TR,"$tRNA_gff" || die "Cannot open the file '$tRNA_gff'.\n";
while(<TR>)
{
	chomp;
	my @sp=split(/\t/,$_);
	if($_ ne "")
	{
		for(my $i=$sp[3];$i<=$sp[4];$i++)
		{
			$pos{$i}=1;
		}
	}
}
close TR;

print "Reading the '$gene_gff' file ...\n";
# CP084950.1	MFAnnot	exon	84	1541	.	+	.	ID=exon:nad4;Parent=mRNA:nad4
open GG,"$gene_gff" || die "Cannot open the file '$gene_gff'.\n";
while(<GG>)
{
	chomp;
	my @sp=split(/\t/,$_);
	if($_=~/\texon\t/)
	{
		for(my $i=$sp[3];$i<=$sp[4];$i++)
		{
			$pos{$i}=1;
		}
	}
}
close GG;

print "Reading the '$DNA_Analyzer_Palindrome' file ...\n";
# Length - Spacer - Mismatches 	Position 	¦¤G(cf) - ¦¤G(lin) 	Palindrome 	Details
# 6-1-0 	54 	2.44 	TTTAAA G TTTAAA 	
open OUT,">$output_file";
open F,"$DNA_Analyzer_Palindrome" || die "Cannot open the file '$DNA_Analyzer_Palindrome'.\n";
while(<F>)
{
	chomp;
	my @sp=split(/\t/,$_);
	if($_=~/^Length/)
	{
		print OUT "$_\n";
	}
	elsif($_ ne "")
	{
		$sp[0]=~s/\s+//g;
		my @rpos=split(/\-/,$sp[0]);
		$sp[2]=~s/\s+//g;
		my $end=$sp[2]+($rpos[0]-1)+$rpos[1]+$rpos[0];
		my $tag=0;
		for(my $i=$sp[2];$i<=$end;$i++)
		{
			if(defined $pos{$i} && $pos{$i}==1)
			{
				$tag=1;
				last;
			}
		}
		if($tag==0)
		{
			print OUT "$_\n";
		}
	}
}
close F;
close OUT;

($timesecond, $timeminute, $timehour, $timedaymonths, $timemonth, $timeyear, $timedayweek, $timedayYear, $timedaylightsavings) = localtime();
$timeyear+=1900;
print "[$timehour\:$timeminute\:$timesecond\, $months[$timemonth] $timedaymonths\, $timeyear] ";
print "End of running the \'$0\' program.\n";

__END__
