#! /usr/bin/perl -w
################################################################################
# fungal_mitogenome_candidate_stop_codon.pl version 1.0
# Copyright (C) Runmao Lin, Tong Liu, Xiaoting Wang, Fanxing Yang, Zhiyin Wang, 2023
# Contact (E-mail): linrunmao@hainanu.cn
# 
# This program is provided under the terms of a personal license to the recipient and 
# may only be used for the recipient's own research at an academic insititution.
# 
# For using this program in a company or for commercial purposes, a commercial license 
# is required.
# 
# The purpose of this program is to identify candidate stop codons for fungal gene.
################################################################################

use strict;
use warnings;
use Getopt::Long;
use FindBin qw($Bin);
use Data::Dumper;

my $gene_gff;
my $gene_id;
my $genome;
my $output;
my $help;
my $version=sprintf("%.1f",1.0);

GetOptions
(
	"gene_gff=s" => \$gene_gff,              # string
	"gene_id=s" => \$gene_id,                # string
	"genome=s" => \$genome,                  # string
	"output=s" => \$output,                  # string
	"help" => \$help                         # flag
);

################################################################################
# usage
############################# usage begin ######################################
my $usage= << "USAGE";

Example: $0 -gene_gff  gene.gff  -gene_id  rps3  -genome genome.fa  -output  rps3_alternative_initiation_codon.txt 
version: 1.0
Options:
	-gene_gff <file>           gene gff, such as 'gene.gff'
	-gene_id <strings>         gene ID, such as 'rps3'
	-genome <file>             genome sequence, such as 'genome.fa'
	-output <file>             output file, such as 'candidate_stop_codon.txt'
	-help                      print help information

USAGE
############################## usage end #######################################

if ($help || !(defined $gene_gff) || !(defined $gene_id) || !(defined $genome) || !(defined $output))
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

my %candidate_stop_codon;
$candidate_stop_codon{"TAA"}="U"; # 
$candidate_stop_codon{"TAG"}="U"; # 

my %seq;
my $genome_id;
open S,"$genome" || die "Cannot open the file '$genome'.\n";
while(<S>)
{
	chomp;
	if($_=~/>/)
	{
		my @sp=split(/\s+/,$_);
		$sp[0]=~s/>//;
		$genome_id=$sp[0];
	}
	elsif($_ ne "")
	{
		$seq{$genome_id}.="$_";
		$seq{$genome_id}=~s/\s+//g;
		$seq{$genome_id}=uc $seq{$genome_id};
	}
}
close S;

open OUT,">$output";
print OUT "GeneID\tGenomeID\tStrand\tBegin\tEnd\tCandidate_Stop_Codon\n";

# T203_mitogenome	MFAnnot	gene	546	1355	.	+	.	ID=gene:cox3;Name=cox3
# T203_mitogenome	MFAnnot	mRNA	546	1355	.	+	.	ID=mRNA:cox3;Name=cox3;Parent=gene:cox3
# T203_mitogenome	MFAnnot	CDS	546	1355	.	+	.	ID=CDS:cox3;Parent=mRNA:cox3
# T203_mitogenome	MFAnnot	exon	546	1355	.	+	.	ID=exon:cox3;Parent=mRNA:cox3

my %gene2begin;
my %gene2end;
my %gene2genome;
my %gene2strand;
my %genepos;
open G,"$gene_gff" || die "Cannot open the file '$gene_gff'.\n";
while(<G>)
{
	chomp;
	my @sp=split(/\t/,$_);
	if($_=~/\tgene\t/)
	{
		my $id=(split(/\;/,$sp[8]))[0];
		$id=~s/ID=gene://;
		$gene2begin{$id}=$sp[3];
		$gene2end{$id}=$sp[4];
		$gene2genome{$id}=$sp[0];
		$gene2strand{$id}=$sp[6];
		for(my $i=$sp[3]-1;$i<$sp[4];$i++)
		{
			$genepos{$sp[0]}{$i}=1;
		}
	}
}
close G;

my @genome_seq=split(//,$seq{$gene2genome{$gene_id}});
if($gene2strand{$gene_id} eq "+")
{
	my $pos=$gene2end{$gene_id}+1000;
	for(my $i=$gene2end{$gene_id};$i<@genome_seq;)
	{
		my $code_seq="$genome_seq[$i]$genome_seq[$i+1]$genome_seq[$i+2]";
		if(defined $genepos{$gene2genome{$gene_id}}{$i} || defined $genepos{$gene2genome{$gene_id}}{$i+1} || defined $genepos{$gene2genome{$gene_id}}{$i+2})
		{
			last;
		}
		if(defined $candidate_stop_codon{$code_seq})
		{
			print OUT "$gene_id\t$gene2genome{$gene_id}\t$gene2strand{$gene_id}\t".($i+1)."\t".($i+3)."\t$candidate_stop_codon{$code_seq}\n";
			last;
		}
		$i+=3;
		last if($i>=$pos);
	}
}
else
{
	my $pos=0;
	if($gene2begin{$gene_id}<1000)
	{
		$pos=1;
	}
	else
	{
		$pos=$gene2begin{$gene_id}-999;
	}
	for(my $i=$gene2begin{$gene_id}-1;$i>$pos;$i--)
	{
		foreach my $gid(keys %gene2genome)
		{
			if($gid ne $gene_id && defined $genepos{$gene2genome{$gid}}{$i})
			{
				$pos=$i+1;
				last;
			}
		}
	}
	my $up_stream=substr($seq{$gene2genome{$gene_id}},$pos-1,$gene2end{$gene_id}-$pos+1);
#	print "***$pos***$gene2end{$gene_id}******\n";
# A: $pos, $pos+1, ..., $gene2begin{$gene_id}, ..., gene2end{$gene_id}
# B: 1, 2, ..., gene2end{$gene_id}-$pos+1
# B(R): gene2end{$gene_id}-$pos+1, gene2end{$gene_id}-$pos, ..., 1
# D: 1(R), 2(R), ..., gene2end{$gene_id}-$pos+1(R)
# n(R)=B:gene2end{$gene_id}-$pos+1-n+1
# n(R)=A:gene2end{$gene_id}-$pos+1-n+1+($pos-1)=gene2end{$gene_id}-n+1
	$up_stream=~tr/ATCG/TAGC/;
	$up_stream=reverse $up_stream;
	my @up_stream_seq=split(//,$up_stream);
#	print "***$up_stream***\n";
	for(my $j=0;$j<@up_stream_seq;)
	{
		my $code_seq="$up_stream_seq[$j]$up_stream_seq[$j+1]$up_stream_seq[$j+2]";
#		my $kk=$gene2end{$gene_id}-$j;
#		print "***$j***$kk***$code_seq***\n";
		if(defined $candidate_stop_codon{$code_seq})
		{
			print OUT "$gene_id\t$gene2genome{$gene_id}\t$gene2strand{$gene_id}\t".($gene2end{$gene_id}-$j-2)."\t".($gene2end{$gene_id}-$j)."\t$candidate_stop_codon{$code_seq}\n";
			last;
		}
		$j+=3;
	}
}
close OUT;

($timesecond, $timeminute, $timehour, $timedaymonths, $timemonth, $timeyear, $timedayweek, $timedayYear, $timedaylightsavings) = localtime();
$timeyear+=1900;
print "[$timehour\:$timeminute\:$timesecond\, $months[$timemonth] $timedaymonths\, $timeyear] ";
print "End of running the \'$0\' program.\n";

__END__
