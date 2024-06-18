#! /usr/bin/perl -w
################################################################################
# select_mitos_MFannot_gff.pl version 1.0
# Copyright (C) Runmao Lin, Tong Liu, Xiaoting Wang, Fanxing Yang, Zhiyin Wang, 2023
# Contact (E-mail): linrunmao@hainanu.cn
#
# This program is provided under the terms of a personal license to the recipient and
# may only be used for the recipient's own research at an academic insititution.
#
# For using this program in a company or for commercial purposes, a commercial license
# is required.
#
# The purpose of this program is to select genes from Mitos and MFanno gff files.
################################################################################

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;

my $select_gene;
my $mitos_update_gff;
my $mfanno_update_gff;
my $output_file;
my $help;

GetOptions
(
	"select_gene=s" => \$select_gene,                                   # string
	"mitos_update_gff=s" => \$mitos_update_gff,                         # string
	"mfanno_update_gff=s" => \$mfanno_update_gff,                       # string
	"output_file=s" => \$output_file,                                   # string
	"help" => \$help                                                    # flag
);

################################################################################
# usage
############################# usage begin ######################################
my $usage= << "USAGE";

Example: perl $0  -select_gene  select_gene.txt  -mitos_update_gff  mitos.update.gene.gff  -mfanno_update_gff  mfanno.update.gene.gff  -output_file  integrated_genes.gff 
version: 1.0
Options:
        -select_gene <file>                    selected gene inforation, such as 'select_gene.txt'
        -mitos_update_gff <file>               mitos gff file, such as 'mitos.update.gene.gff'
        -mfanno_update_gff <file>              mfanno gff file, such as 'mfanno.update.gene.gff'
        -output_file <file>                    output file, such as 'infor.txt'
        -help                                  print help information

USAGE
############################## usage end #######################################

if ($help || !(defined $select_gene) || !(defined $mitos_update_gff) || !(defined $mfanno_update_gff) || !(defined $output_file))
{
	print $usage;
	exit;
}

my @months = qw(Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec);
my ($timesecond, $timeminute, $timehour, $timedaymonths, $timemonth, $timeyear, $timedayweek, $timedayYear, $timedaylightsavings) = localtime();
$timeyear+=1900;

print "[$timehour\:$timeminute\:$timesecond\, $months[$timemonth] $timedaymonths\, $timeyear] ";
print "Start to run the \'$0\' program ...\n";

############################## select_gene.txt #################################
# gff_file	GeneID	NewID
# mfanno.update.gene.gff	atp6	atp6
# mfanno.update.gene.gff	cox2_2	cox2
# result.update.gene.gff	cob_0	cob
################################################################################

################################### gene.gff ###################################
# MW525445.1	mitos	gene	2090	87373	.	+	.	ID=gene:giy_11;Name=giy_11
# MW525445.1	mitos	mRNA	2090	87373	.	+	.	ID=mRNA:giy_11;Name=giy_11;Parent=gene:giy_11
# MW525445.1	mitos	CDS	2090	2431	.	+	.	ID=CDS:giy_11;Parent=mRNA:giy_11
# MW525445.1	mitos	exon	2090	2431	.	+	.	ID=exon:giy_11;Parent=mRNA:giy_11
# MW525445.1	mitos	CDS	87167	87373	.	+	.	ID=CDS:giy_11;Parent=mRNA:giy_11
# MW525445.1	mitos	exon	87167	87373	.	+	.	ID=exon:giy_11;Parent=mRNA:giy_11
################################################################################

print "Reading the '$mitos_update_gff' file ...\n";
my $mitos_gene=reading_gff($mitos_update_gff);

print "Reading the '$mfanno_update_gff' file ...\n";
my $mfanno_gene=reading_gff($mfanno_update_gff);

my %select_gene;
print "Reading the '$select_gene' file ...\n";
open F,"$select_gene" || die "Cannot open the file '$select_gene'.\n";
while(<F>)
{
	chomp;
	my @sp=split(/\t/,$_);
	if($_ ne "" && !($_=~/^gff_file/) && $sp[0]=~/mfanno/i)
	{
		if(!(defined $mfanno_gene->{$sp[1]}))
		{
			print "Error! No '$sp[1]' in '$sp[0]' or '$mfanno_update_gff'. Please check!\n";
			exit;
		}
		if($sp[1] ne $sp[2])
		{
			my $gene_ID_gff=$mfanno_gene->{$sp[1]};
			$gene_ID_gff=~s/ID=gene:$sp[1]/ID=gene:$sp[2]/g;
			$gene_ID_gff=~s/Name=$sp[1]/Name=$sp[2]/g;
			$gene_ID_gff=~s/ID=mRNA:$sp[1]/ID=mRNA:$sp[2]/g;
			$gene_ID_gff=~s/Parent=gene:$sp[1]/Parent=gene:$sp[2]/g;
			$gene_ID_gff=~s/ID=CDS:$sp[1]/ID=CDS:$sp[2]/g;
			$gene_ID_gff=~s/Parent=mRNA:$sp[1]/Parent=mRNA:$sp[2]/g;
			$gene_ID_gff=~s/ID=exon:$sp[1]/ID=exon:$sp[2]/g;
			my @gsp=split(/\n/,$gene_ID_gff);
			my @ssp=split(/\t/,$gsp[0]);
			$select_gene{$ssp[0]}{$ssp[3]}=$gene_ID_gff;
		}
		else
		{
			my @gsp=split(/\n/,$mfanno_gene->{$sp[1]});
			my @ssp=split(/\t/,$gsp[0]);
			$select_gene{$ssp[0]}{$ssp[3]}=$mfanno_gene->{$sp[1]};
		}
	}
	elsif($_ ne "" && !($_=~/^gff_file/) && ($sp[0]=~/mitos/i || $sp[0]=~/result.update.gene.gff/))
	{
		if(!(defined $mitos_gene->{$sp[1]}))
		{
			print "Error! No '$sp[1]' in '$sp[0]' or '$mitos_update_gff'. Please check!\n";
			exit;
		}
		if($sp[1] ne $sp[2])
		{
			my $gene_ID_gff=$mitos_gene->{$sp[1]};
			$gene_ID_gff=~s/ID=gene:$sp[1]/ID=gene:$sp[2]/g;
			$gene_ID_gff=~s/Name=$sp[1]/Name=$sp[2]/g;
			$gene_ID_gff=~s/ID=mRNA:$sp[1]/ID=mRNA:$sp[2]/g;
			$gene_ID_gff=~s/Parent=gene:$sp[1]/Parent=gene:$sp[2]/g;
			$gene_ID_gff=~s/ID=CDS:$sp[1]/ID=CDS:$sp[2]/g;
			$gene_ID_gff=~s/Parent=mRNA:$sp[1]/Parent=mRNA:$sp[2]/g;
			$gene_ID_gff=~s/ID=exon:$sp[1]/ID=exon:$sp[2]/g;
			my @gsp=split(/\n/,$gene_ID_gff);
			my @ssp=split(/\t/,$gsp[0]);
			$select_gene{$ssp[0]}{$ssp[3]}=$gene_ID_gff;
		}
		else
		{
			my @gsp=split(/\n/,$mitos_gene->{$sp[1]});
			my @ssp=split(/\t/,$gsp[0]);
			$select_gene{$ssp[0]}{$ssp[3]}=$mitos_gene->{$sp[1]};
		}
	}
}
close F;

open OUT,">$output_file";
foreach my $sid(sort keys %select_gene)
{
	foreach my $gpos(sort {$a<=>$b} keys %{$select_gene{$sid}})
	{
		print OUT $select_gene{$sid}{$gpos};
	}
}
close OUT;


sub reading_gff
{
	my $gff_file=shift;
	my $gene_infor;
	open G,"$gff_file" || die "Cannot open the file '$gff_file'.\n";
	while(<G>)
	{
		chomp;
		my @sp=split(/\t/,$_);
		if($_=~/\tgene\t/)
		{
			my $id=(split(/\;/,$sp[8]))[-1];
			$id=~s/Name=//;
			$gene_infor->{$id}="$_\n";
		}
		elsif($_=~/\tmRNA\t/)
		{
			my $id=(split(/\;/,$sp[8]))[-1];
			$id=~s/Parent=gene://;
			$gene_infor->{$id}.="$_\n";
		}
		elsif($_=~/\tCDS\t/)
		{
			my $id=(split(/\;/,$sp[8]))[-1];
			$id=~s/Parent=mRNA://;
			$gene_infor->{$id}.="$_\n";
		}
		elsif($_=~/\texon\t/)
		{
			my $id=(split(/\;/,$sp[8]))[-1];
			$id=~s/Parent=mRNA://;
			$gene_infor->{$id}.="$_\n";
		}
	}
	close G;
	return $gene_infor;
}

($timesecond, $timeminute, $timehour, $timedaymonths, $timemonth, $timeyear, $timedayweek, $timedayYear, $timedaylightsavings) = localtime();
$timeyear+=1900;
print "[$timehour\:$timeminute\:$timesecond\, $months[$timemonth] $timedaymonths\, $timeyear] ";
print "End of running the \'$0\' program.\n";

__END__
