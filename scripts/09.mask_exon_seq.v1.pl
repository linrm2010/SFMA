#! /usr/bin/perl -w
################################################################################
# mask_exon_seq.pl version 1.0
# Copyright (C) Runmao Lin, Tong Liu, Xiaoting Wang, Fanxing Yang, Zhiyin Wang, 2023
# Contact (E-mail): linrunmao@hainanu.cn
#
# This program is provided under the terms of a personal license to the recipient and
# may only be used for the recipient's own research at an academic insititution.
#
# For using this program in a company or for commercial purposes, a commercial license
# is required.
#
# The purpose of this program is to mask exon sequences.
################################################################################

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;

my $genome;
my $gene_gff;
my $tRNA_gff;
my $rRNA_gff;
my $output_file;
my $help;

GetOptions
(
	"genome=s" => \$genome,                         # string
	"gene_gff=s" => \$gene_gff,                     # string
	"tRNA_gff=s" => \$tRNA_gff,                     # string
	"rRNA_gff=s" => \$rRNA_gff,                     # string
	"output_file=s" => \$output_file,               # string
	"help" => \$help                                # flag
);

################################################################################
# usage
############################# usage begin ######################################
my $usage= << "USAGE";

Example: perl $0 -genome  genome.fa  -gene_gff  integrate.gene.gff  -tRNA_gff  tRNA_result.gff  -rRNA_gff  rRNA.gff  -output_file  genome.mask_genes.fa 
version: 1.0
Options:
        -genome <file>               genome sequence file, such as 'genome.fa'
        -gene_gff <file>             gene gff file, such as 'integrate.gene.gff'
        -tRNA_gff <file>             tRNA gff file, such as 'tRNA_result.gff'
        -rRNA_gff <file>             rRNA gff file, such as 'rRNA.gff'
        -output_file <file>          output file, such as 'genome.mask_genes.fa'
        -help                        print help information

USAGE
############################## usage end #######################################

if ($help || !(defined $genome) || !(defined $gene_gff) || !(defined $tRNA_gff) || !(defined $rRNA_gff) || !(defined $output_file))
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
open S,"$genome" || die "Cannot open the file '$genome'.\n";
$/=">";
<S>;
$/="\n";
while(<S>)
{
	chomp;
	my $id=(split(/\s+/,$_))[0];
	$/=">";
	$seq{$id}=<S>;
	$seq{$id}=~s/>//g;
	$seq{$id}=~s/\s+//g;
	$/="\n";
}
$/="\n";
close S;

################################################################################
# Reading 'gene_gff'
################################################################################

print "Reading '$gene_gff' file ...\n";

my %gene;
open G,"$gene_gff" || die "Cannot open the file '$gene_gff'.\n";
while(<G>)
{
	chomp;
	my @sp=split(/\t/,$_);
	if($_=~/\texon\t/)
	{
		$gene{"$sp[0]\t$sp[3]\t$sp[4]"}=1;
	}
}
close G;

################################################################################
# Reading 'tRNA_gff'
################################################################################

print "Reading '$tRNA_gff' file ...\n";

open T,"$tRNA_gff" || die "Cannot open the file '$tRNA_gff'.\n";
while(<T>)
{
	chomp;
	my @sp=split(/\t/,$_);
	if($_=~/\ttRNA\t/)
	{
		$gene{"$sp[0]\t$sp[3]\t$sp[4]"}=1;
	}
}
close T;

################################################################################
# Reading 'rRNA_gff'
################################################################################

print "Reading '$rRNA_gff' file ...\n";

open R,"$rRNA_gff" || die "Cannot open the file '$rRNA_gff'.\n";
while(<R>)
{
	chomp;
	my @sp=split(/\t/,$_);
	if($_=~/\trRNA\t/)
	{
		$gene{"$sp[0]\t$sp[3]\t$sp[4]"}=1;
	}
}
close R;

open OUT,">$output_file";
foreach my $id(keys %seq)
{
	foreach my $pid(keys %gene)
	{
		my @sp=split(//,$seq{$id});
		my $seq_update="";
		my @psp=split(/\t/,$pid);
		if($psp[0] eq $id)
		{
			for(my $i=0;$i<($psp[1]-1);$i++)
			{
				$seq_update.=$sp[$i];
			}
			for(my $i=$psp[1]-1;$i<$psp[2];$i++)
			{
				$seq_update.="N";
			}
			for(my $i=$psp[2];$i<@sp;$i++)
			{
				$seq_update.=$sp[$i];
			}
			$seq{$id}=$seq_update;
		}
	}
	my $seq_reset="";
	my @sp=split(//,$seq{$id});
	for(my $i=0;$i<@sp;$i++)
	{
		$seq_reset.="$sp[$i]";
		$seq_reset.="\n" if(($i+1)%100==0);
	}
	$seq_reset=~s/\n$//;
	print OUT ">$id\n$seq_reset\n";
}
close OUT;

($timesecond, $timeminute, $timehour, $timedaymonths, $timemonth, $timeyear, $timedayweek, $timedayYear, $timedaylightsavings) = localtime();
$timeyear+=1900;
print "[$timehour\:$timeminute\:$timesecond\, $months[$timemonth] $timedaymonths\, $timeyear] ";
print "End of running the \'$0\' program.\n";

__END__
