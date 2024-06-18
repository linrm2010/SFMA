#! /usr/bin/perl -w
################################################################################
# comparison_mitos_MFannot.pl version 1.0
# Copyright (C) Runmao Lin, Tong Liu, Xiaoting Wang, Fanxing Yang, Zhiyin Wang, 2023
# Contact (E-mail): linrunmao@hainanu.cn
#
# This program is provided under the terms of a personal license to the recipient and
# may only be used for the recipient's own research at an academic insititution.
#
# For using this program in a company or for commercial purposes, a commercial license
# is required.
#
# The purpose of this program is to compare predicted results from Mitos and MFanno.
################################################################################

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;

my $mitos_update_gff;
my $mfanno_update_gff;
my $output_file;
my $help;

GetOptions
(
	"mitos_update_gff=s" => \$mitos_update_gff,                         # string
	"mfanno_update_gff=s" => \$mfanno_update_gff,                       # string
	"output_file=s" => \$output_file,                                   # string
	"help" => \$help                                                    # flag
);

################################################################################
# usage
############################# usage begin ######################################
my $usage= << "USAGE";

Example: perl $0 -mitos_update_gff  mitos.update.gene.gff  -mfanno_update_gff  mfanno.update.gene.gff  -output_file  infor.txt 
version: 1.0
Options:
        -mitos_update_gff <file>               mitos gff file, such as 'mitos.update.gene.gff'
        -mfanno_update_gff <file>              mfanno gff file, such as 'mfanno.update.gene.gff'
        -output_file <file>                    output file, such as 'infor.txt'
        -help                                  print help information

USAGE
############################## usage end #######################################

if ($help || !(defined $mitos_update_gff) || !(defined $mfanno_update_gff) || !(defined $output_file))
{
	print $usage;
	exit;
}

my @months = qw(Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec);
my ($timesecond, $timeminute, $timehour, $timedaymonths, $timemonth, $timeyear, $timedayweek, $timedayYear, $timedaylightsavings) = localtime();
$timeyear+=1900;

print "[$timehour\:$timeminute\:$timesecond\, $months[$timemonth] $timedaymonths\, $timeyear] ";
print "Start to run the \'$0\' program ...\n";

# mitos
# MW525445.1	mitos	gene	5091	8774	.	+	.	ID=gene:cox3;Name=cox3
# MW525445.1	mitos	mRNA	5091	8774	.	+	.	ID=mRNA:cox3;Name=cox3;Parent=gene:cox3
# MW525445.1	mitos	CDS	5091	5729	.	+	.	ID=CDS:cox3;Parent=mRNA:cox3
# MW525445.1	mitos	exon	5091	5729	.	+	.	ID=exon:cox3;Parent=mRNA:cox3
# MW525445.1	mitos	CDS	8604	8774	.	+	.	ID=CDS:cox3;Parent=mRNA:cox3
# MW525445.1	mitos	exon	8604	8774	.	+	.	ID=exon:cox3;Parent=mRNA:cox3

# mfanno
# MW525445.1	MFAnnot	gene	5091	8774	.	+	.	ID=gene:cox3;Name=cox3
# MW525445.1	MFAnnot	mRNA	5091	8774	.	+	.	ID=mRNA:cox3;Name=cox3;Parent=gene:cox3
# MW525445.1	MFAnnot	CDS	5091	5730	.	+	.	ID=CDS:cox3;Parent=mRNA:cox3
# MW525445.1	MFAnnot	CDS	8605	8774	.	+	.	ID=CDS:cox3;Parent=mRNA:cox3
# MW525445.1	MFAnnot	exon	5091	5730	.	+	.	ID=exon:cox3;Parent=mRNA:cox3
# MW525445.1	MFAnnot	exon	8605	8774	.	+	.	ID=exon:cox3;Parent=mRNA:cox3

my $mitos_gene;
my $mitos_cds;
my $mitos_length;
my $mfanno_gene;
my $mfanno_cds;
my $mfanno_length;

################################################################################
# Reading 'mitos_update_gff'
################################################################################

print "Reading '$mitos_update_gff' file ...\n";

($mitos_gene,$mitos_cds,$mitos_length)=reading_gff($mitos_update_gff);

################################################################################
# Reading 'mfanno_update_gff'
################################################################################

print "Reading '$mfanno_update_gff' file ...\n";

($mfanno_gene,$mfanno_cds,$mfanno_length)=reading_gff($mfanno_update_gff);

my $mitos_mfanno;
my $mitos_only;
my $mfanno_only;
my %check_mfanno_gene;
foreach my $mid(sort keys %{$mitos_gene})
{
	my $mid_unique=1;
	foreach my $fid(sort keys %{$mfanno_gene})
	{
		if($mitos_gene->{$mid} eq $mfanno_gene->{$fid} && $mitos_cds->{$mid} eq $mfanno_cds->{$fid})
		{
			$mitos_mfanno.="$mid\, $fid\n";
			$mid_unique=0;
			$check_mfanno_gene{$fid}=1;
			last;
		}
	}
	if($mid_unique==1 && ($mitos_length->{$mid}/3)>39)
	{
		$mitos_only.="$mid\n";
	}
	elsif($mid_unique==1 && ($mitos_length->{$mid}/3)<=39)
	{
		$mitos_only.="$mid \| ignore; aa length < 40\n";
	}
}
foreach my $fid(sort keys %{$mfanno_gene})
{
	if(!(defined $check_mfanno_gene{$fid}) && ($mfanno_length->{$fid}/3)>39)
	{
		$mfanno_only.="$fid\n";
	}
	elsif(!(defined $check_mfanno_gene{$fid}) && ($mfanno_length->{$fid}/3)<=39)
	{
		$mfanno_only.="$fid \| ignore; aa length < 40\n";
	}
}

open OUT,">$output_file";
print OUT "# mitos and mfannot\n";
if(defined $mitos_mfanno && $mitos_mfanno ne "")
{
	print OUT $mitos_mfanno;
}
print OUT "\n";
print OUT "# mitos unique\n";
if(defined $mitos_only && $mitos_only ne "")
{
	print OUT $mitos_only;
}
print OUT "\n";
print OUT "# mfanno unique\n";
if(defined $mfanno_only && $mfanno_only ne "")
{
	print OUT $mfanno_only;
}
close OUT;


sub reading_gff
{
	my $gff_file=shift;
	my $gene_infor;
	my $cds_infor;
	my $cds_length;
	open G,"$gff_file" || die "Cannot open the file '$gff_file'.\n";
	while(<G>)
	{
		chomp;
		my @sp=split(/\t/,$_);
		if($_=~/\tgene\t/)
		{
			my $id=(split(/\;/,$sp[8]))[-1];
			$id=~s/Name=//;
			$gene_infor->{$id}="$sp[0]\t$sp[3]\t$sp[4]";
		}
		elsif($_=~/\tCDS\t/)
		{
			my $id=(split(/\;/,$sp[8]))[0];
			$id=~s/ID=CDS\://;
			$cds_infor->{$id}.="$sp[0]\t$sp[3]\t$sp[4]\;";
			$cds_length->{$id}+=$sp[4]-$sp[3]+1;
		}
	}
	close G;
	return ($gene_infor,$cds_infor,$cds_length);
}

($timesecond, $timeminute, $timehour, $timedaymonths, $timemonth, $timeyear, $timedayweek, $timedayYear, $timedaylightsavings) = localtime();
$timeyear+=1900;
print "[$timehour\:$timeminute\:$timesecond\, $months[$timemonth] $timedaymonths\, $timeyear] ";
print "End of running the \'$0\' program.\n";

__END__
