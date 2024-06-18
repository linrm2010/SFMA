#! /usr/bin/perl -w
################################################################################
# fungal_mitogenome_gff2cds_check_complete.pl version 1.0
# Copyright (C) Runmao Lin, Tong Liu, Xiaoting Wang, Fanxing Yang, Zhiyin Wang, 2023
# Contact (E-mail): linrunmao@hainanu.cn
#
# This program is provided under the terms of a personal license to the recipient and
# may only be used for the recipient's own research at an academic insititution.
#
# For using this program in a company or for commercial purposes, a commercial license
# is required.
#
# The purpose of this program is to extract cds sequences.
# The output results will be placed on current direction.
################################################################################

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;

my $genome;
my $gff_file;
my $output_file;
my $help;

GetOptions
(
	"genome=s" => \$genome,                                             # string
	"gff_file=s" => \$gff_file,                                         # string
	"output_file=s" => \$output_file,                                   # string
	"help" => \$help                                                    # flag
);

################################################################################
# usage
############################# usage begin ######################################
my $usage= << "USAGE";

Example: perl $0 -genome  genome.fa  -gff_file  gene.gff  -output_file  gene.cds 
version: 1.0
Options:
        -genome <file>                         genome sequences, such as 'genome.fa'
        -gff_file <file>                       the gene gff file, such as 'gene.gff'
        -output_file <file>                    output file, such as 'gene.cds'
        -help                                  print help information

USAGE
############################## usage end #######################################

if ($help || !(defined $genome) || !(defined $gff_file) || !(defined $output_file))
{
	print $usage;
	exit;
}

my @months = qw(Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec);
my ($timesecond, $timeminute, $timehour, $timedaymonths, $timemonth, $timeyear, $timedayweek, $timedayYear, $timedaylightsavings) = localtime();
$timeyear+=1900;

print "[$timehour\:$timeminute\:$timesecond\, $months[$timemonth] $timedaymonths\, $timeyear] ";
print "Start to run the \'$0\' program ...\n";

## case: if there are more 'three_prime_UTR' or 'five_prime_UTR'
## contig24	EuGene	gene	155078	156178	.	-	.	ID=gene:AG1IC_00353;Name=AG1IC_00353;length=1101
## contig24	EuGene	mRNA	155078	156178	.	-	.	ID=mRNA:AG1IC_00353;Name=AG1IC_00353;Parent=gene:AG1IC_00353;nb_exon=4;length=875
## contig24	EuGene	exon	155078	155111	.	-	.	ID=exon:AG1IC_00353.0;Parent=mRNA:AG1IC_00353;Ontology_term=SO:0000205
## contig24	EuGene	exon	155163	155307	.	-	.	ID=exon:AG1IC_00353.4;Parent=mRNA:AG1IC_00353;Ontology_term=SO:0000202
## contig24	EuGene	exon	155356	155728	.	-	1	ID=exon:AG1IC_00353.3;Parent=mRNA:AG1IC_00353;Ontology_term=SO:0000004
## contig24	EuGene	exon	155803	155920	.	-	2	ID=exon:AG1IC_00353.2;Parent=mRNA:AG1IC_00353;Ontology_term=SO:0000004
## contig24	EuGene	exon	155974	156178	.	-	0	ID=exon:AG1IC_00353.1;Parent=mRNA:AG1IC_00353;Ontology_term=SO:0000200
## contig24	EuGene	three_prime_UTR	155078	155111	.	-	.	ID=three_prime_UTR:AG1IC_00353.0;Parent=mRNA:AG1IC_00353;Ontology_term=SO:0000205
## contig24	EuGene	three_prime_UTR	155163	155250	.	-	.	ID=three_prime_UTR:AG1IC_00353.2;Parent=mRNA:AG1IC_00353;Ontology_term=SO:0000205
## contig24	EuGene	CDS	155251	155307	.	-	0	ID=CDS:AG1IC_00353.4;Parent=mRNA:AG1IC_00353;Ontology_term=SO:0000197
## contig24	EuGene	CDS	155356	155728	.	-	1	ID=CDS:AG1IC_00353.3;Parent=mRNA:AG1IC_00353;Ontology_term=SO:0000004
## contig24	EuGene	CDS	155803	155920	.	-	2	ID=CDS:AG1IC_00353.2;Parent=mRNA:AG1IC_00353;Ontology_term=SO:0000004
## contig24	EuGene	CDS	155974	156148	.	-	0	ID=CDS:AG1IC_00353.1;Parent=mRNA:AG1IC_00353;Ontology_term=SO:0000196
## contig24	EuGene	five_prime_UTR	156149	156178	.	-	.	ID=five_prime_UTR:AG1IC_00353.10;Parent=mRNA:AG1IC_00353;Ontology_term=SO:0000204

my %seq;
open S,"$genome" || die "Cannot open the file '$genome'.\n";
$/=">";<S>;$/="\n";
while(<S>)
{
	my $id=$1 if(/(\S+)/);
	$/=">";
	$seq{$id}=<S>;
	chomp $seq{$id};
	$seq{$id}=~s/\s+//g;
	$seq{$id}=uc $seq{$id};
	$/="\n";
}
close S;

my %cds;
my %strand;
my %gene;
open IN,"$gff_file" || die "Cannot open the file '$gff_file'.\n";
while(<IN>)
{
	chomp;
	my @sp=split(/\t/,$_);
	if($_=~/\tgene\t/)
	{
		my $geneid=(split(/\;/,$sp[8]))[0];
		$geneid=~s/ID=gene\://;
		$gene{$geneid}="$sp[0]\:$sp[3]\:$sp[4]\:$sp[6]";
		$strand{$geneid}=$sp[6];
	}
	if($_=~/\tCDS\t/)
	{
		my $geneid=(split(/\;/,$sp[8]))[0];
		$geneid=~s/ID=CDS\://;
		$geneid=~s/\.\d+$//;
		$cds{$geneid}.=substr($seq{$sp[0]},$sp[3]-1,$sp[4]-$sp[3]+1);
	}
}
close IN;

open OUT,">$output_file";
open OUTERR,">$output_file\.err";
foreach my $id(sort keys %cds)
{
	if($strand{$id} eq "-")
	{
		$cds{$id}=~tr/ATCG/TAGC/;
		$cds{$id}=reverse $cds{$id};
	}
	print OUTERR "$id\n" if($cds{$id} eq "");
	if(!($cds{$id}=~/^ATG/))
	{
		print OUTERR "$id Not complete, without begin 'ATG' in '$id'.\n";
	}
	if( !($cds{$id}=~/TAA$/ || $cds{$id}=~/TAG$/))
	{
		print OUTERR "$id Not complete, without end 'TAA/TAG' in '$id'.\n";
	}
	my @csp=split(//,$cds{$id});
	$cds{$id}="";
	for(my $i=0;$i<@csp;$i++)
	{
		$cds{$id}.=$csp[$i];
		$cds{$id}.="\n" if(($i+1)%60==0);
	}
	$cds{$id}=~s/\n$//;
	print OUT "\>$id   $gene{$id}\n$cds{$id}\n";
}
close OUT;
close OUTERR;

($timesecond, $timeminute, $timehour, $timedaymonths, $timemonth, $timeyear, $timedayweek, $timedayYear, $timedaylightsavings) = localtime();
$timeyear+=1900;
print "[$timehour\:$timeminute\:$timesecond\, $months[$timemonth] $timedaymonths\, $timeyear] ";
print "End of running the \'$0\' program.\n";

__END__
