#! /usr/bin/perl -w
################################################################################
# fungal_mitogenome_alternative_initiation_codon.pl version 1.0
# Copyright (C) Runmao Lin, Tong Liu, Xiaoting Wang, Fanxing Yang, Zhiyin Wang, 2023
# Contact (E-mail): linrunmao@hainanu.cn
# 
# This program is provided under the terms of a personal license to the recipient and 
# may only be used for the recipient's own research at an academic insititution.
# 
# For using this program in a company or for commercial purposes, a commercial license 
# is required.
# 
# The purpose of this program is to identify alternative initiation codons in fungal gene.
################################################################################

use strict;
use warnings;
use Getopt::Long;
use FindBin qw($Bin);
use Data::Dumper;

my $gene_cds;
my $gene_id;
my $gene_strand;
my $output;
my $help;
my $version=sprintf("%.1f",1.0);

############################## gene cds format #################################
#>atp8 
#ATGCCACAATTAACACCTTTTTATTATATGAACGAAATAGTATTTGCTTTTGCAATAATA
#GTAATACTATTATTCATACTATCTAAATACATACTACCAAGAATAGTACGTTTATTTTTA
#TCACGTATGTTTATTAATAAAATATAA
################################################################################

GetOptions
(
	"gene_cds=s" => \$gene_cds,              # string
	"gene_id=s" => \$gene_id,                # string
	"gene_strand=s" => \$gene_strand,        # string
	"output=s" => \$output,                  # string
	"help" => \$help                         # flag
);

################################################################################
# usage
############################# usage begin ######################################
my $usage= << "USAGE";

Example: $0 -gene_cds  gene.cds  -gene_id  rps3  -gene_strand +  -output  rps3_alternative_initiation_codon.txt 
version: 1.0
Options:
	-gene_cds <file>           one file contains gene nucleotide sequences, such as 'gene.cds'
	-gene_id <strings>         gene ID, such as 'rps3'
	-gene_strand <strings>     strand of coding genes, such as '+'
	-output <file>             output file, such as 'rps3_alternative_initiation_codon.txt'
	-help                      print help information

USAGE
############################## usage end #######################################

if ($help || !(defined $gene_cds) || !(defined $gene_id) || !(defined $gene_strand) || !(defined $output))
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
# Reading 'gene_cds'
################################################################################

print "Reading '$gene_cds' file ...\n";

my %alternative_IC;
$alternative_IC{"TTA"}="M"; # Trypanosoma
$alternative_IC{"TTG"}="M"; # Trypanosoma
$alternative_IC{"CTG"}="M"; # Trypanosoma
$alternative_IC{"ATT"}="M"; # Leishmania, Tertrahymena, Paramecium
$alternative_IC{"ATA"}="M"; # Leishmania, Tertrahymena, Paramecium
$alternative_IC{"ATG"}="M"; # Tertrahymena, Paramecium
$alternative_IC{"ATC"}="M"; # Paramecium
$alternative_IC{"GTG"}="M"; # Paramecium

my %seq;
my $geneid;
open S,"$gene_cds" || die "Cannot open the file '$gene_cds'.\n";
while(<S>)
{
	chomp;
	if($_=~/>/)
	{
		my @sp=split(/\s+/,$_);
		$sp[0]=~s/>//;
		$geneid=$sp[0];
	}
	elsif($_ ne "")
	{
		$seq{$geneid}.="$_";
		$seq{$geneid}=~s/\s+//g;
		$seq{$geneid}=uc $seq{$geneid};
	}
}
close S;

open OUT,">$output";
print OUT "GeneLength\tCodon_seq\tAlternative_Initiation_Codon\tStrand\tPosition\n";
if(!(defined $seq{$gene_id}))
{
	print "No gene '$gene_id' in '$gene_cds'.\n";
	exit;
}
if((length $seq{$gene_id})%3!=0)
{
	print "Codon error for '$gene_id'?! Length of '$gene_id' in '$gene_cds' was not 3x.\n";
	exit;
}

my $pos=0;
my @gene_seq=split(//,$seq{$gene_id});
for(my $i=0;$i<@gene_seq;)
{
	my $code_seq="$gene_seq[$i]$gene_seq[$i+1]$gene_seq[$i+2]";
	my $strand_infor="";
	if($gene_strand eq "+")
	{
		$strand_infor="$gene_strand\t\+$pos";
	}
	else
	{
		$strand_infor="$gene_strand\t\-$pos";
	}
	if(defined $alternative_IC{$code_seq} && $alternative_IC{$code_seq} eq "M")
	{
		print OUT (scalar @gene_seq)."\t$code_seq\t$alternative_IC{$code_seq}\t$strand_infor\n";
	}
	$pos+=3;
	$i+=3;
}
close OUT;

($timesecond, $timeminute, $timehour, $timedaymonths, $timemonth, $timeyear, $timedayweek, $timedayYear, $timedaylightsavings) = localtime();
$timeyear+=1900;
print "[$timehour\:$timeminute\:$timesecond\, $months[$timemonth] $timedaymonths\, $timeyear] ";
print "End of running the \'$0\' program.\n";

__END__
