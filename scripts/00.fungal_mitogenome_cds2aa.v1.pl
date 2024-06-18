#! /usr/bin/perl -w
################################################################################
# fungal_mitogenome_cds2aa.pl version 1.0
# Copyright (C) Runmao Lin, Tong Liu, Xiaoting Wang, Fanxing Yang, Zhiyin Wang, 2023
# Contact (E-mail): linrunmao@hainanu.cn
#
# This program is provided under the terms of a personal license to the recipient and
# may only be used for the recipient's own research at an academic insititution.
#
# For using this program in a company or for commercial purposes, a commercial license
# is required.
#
# The purpose of this program is to translate nucleotide sequences to amino acids.
################################################################################

use strict;
use warnings;
use Getopt::Long;
use FindBin qw($Bin);
use Data::Dumper;

my $cds_file;
my $aa_file;
my $help;

GetOptions
(
	"cds_file=s" => \$cds_file,                                         # string
	"aa_file=s" => \$aa_file,                                           # string
	"help" => \$help                                                    # flag
);

################################################################################
# usage
############################# usage begin ######################################
my $usage= << "USAGE";

Example: perl $0 -cds_file  cds.fa  -aa_file  pep.fa 
version: 1.0
Options:
        -cds_file <file>                       cds sequences, such as 'cds.fa'
        -aa_file <file>                        output file for amino acids, such as 'pep.fa'
        -help                                  print help information

USAGE
############################## usage end #######################################

if ($help || !(defined $cds_file) || !(defined $aa_file))
{
	print $usage;
	exit;
}

my @months = qw(Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec);
my ($timesecond, $timeminute, $timehour, $timedaymonths, $timemonth, $timeyear, $timedayweek, $timedayYear, $timedaylightsavings) = localtime();
$timeyear+=1900;

print "[$timehour\:$timeminute\:$timesecond\, $months[$timemonth] $timedaymonths\, $timeyear] ";
print "Start to run the \'$0\' program ...\n";

# The Genetic Codes
# https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
# 4. The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code (transl_table=4)
#  AAs  = FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
#Starts = --MM------**-------M------------MMMM---------------M------------
#Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
#Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
#Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG

my %codon;
my $aa=   "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
my $base1="TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG";
my $base2="TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG";
my $base3="TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG";
my @spaa=split(//,$aa);
my @spbase1=split(//,$base1);
my @spbase2=split(//,$base2);
my @spbase3=split(//,$base3);
for(my $i=0;$i<@spbase1;$i++)
{
	$spaa[$i]="U" if($spaa[$i] eq "*");
	$codon{"$spbase1[$i]$spbase2[$i]$spbase3[$i]"}=$spaa[$i];
}

open OUT,">$aa_file";

################################################################################
# Reading 'cds_file'
################################################################################

print "Reading '$cds_file' file ...\n";

open S,"$cds_file" || die "Cannot open the file '$cds_file'.\n";
$/=">";
<S>;
$/="\n";
while(<S>)
{
	my $id=(split(/\s+/,$_))[0];
	$id=~s/\>//g;
	$/=">";
	my $seq=<S>;
	$seq=~s/>//g;
	$seq=~s/\s+//g;
	$seq=uc $seq;
	$/="\n";
	my $aaseq="";
	my @sp=split(//,$seq);
	for(my $i=0;$i<@sp-2;)
	{
		if(!(defined $codon{"$sp[$i]$sp[$i+1]$sp[$i+2]"}))
		{
			$aaseq.="X";
		}
		elsif(!($codon{"$sp[$i]$sp[$i+1]$sp[$i+2]"} eq "U" && ($i+3)==(scalar @sp)))
		{
			$aaseq.=$codon{"$sp[$i]$sp[$i+1]$sp[$i+2]"};
		}
		elsif($codon{"$sp[$i]$sp[$i+1]$sp[$i+2]"} eq "U" && ($i+3)==(scalar @sp))
		{
			$aaseq.="*";
		}
		$aaseq.="\n" if(($i+3)%50==0);
		$i+=3;
	}
	$aaseq=~s/\s+$//;
	print OUT ">$id\n$aaseq\n";
}
$/="\n";
close S;
close OUT;

($timesecond, $timeminute, $timehour, $timedaymonths, $timemonth, $timeyear, $timedayweek, $timedayYear, $timedaylightsavings) = localtime();
$timeyear+=1900;
print "[$timehour\:$timeminute\:$timesecond\, $months[$timemonth] $timedaymonths\, $timeyear] ";
print "End of running the \'$0\' program.\n";

__END__
