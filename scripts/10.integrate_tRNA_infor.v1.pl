#! /usr/bin/perl -w
################################################################################
# integrate_tRNA_infor.pl version 1.0
# Copyright (C) Runmao Lin, Tong Liu, Xiaoting Wang, Fanxing Yang, Zhiyin Wang, 2023
# Contact (E-mail): linrunmao@hainanu.cn
#
# This program is provided under the terms of a personal license to the recipient and
# may only be used for the recipient's own research at an academic insititution.
#
# For using this program in a company or for commercial purposes, a commercial license
# is required.
#
# The purpose of this program is to integrate tRNA genes.
################################################################################

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;

my $mitos_tRNA_gff;
my $mfannot_tRNA_gff;
my $tRNAscan_out;
my $output_prefix;
my $help;

GetOptions
(
	"mitos_tRNA_gff=s" => \$mitos_tRNA_gff,                             # string
	"mfannot_tRNA_gff=s" => \$mfannot_tRNA_gff,                         # string
	"tRNAscan_out=s" => \$tRNAscan_out,                                 # string
	"output_prefix=s" => \$output_prefix,                               # string
	"help" => \$help                                                    # flag
);

################################################################################
# usage
############################# usage begin ######################################
my $usage= << "USAGE";

Example: perl $0  -mitos_tRNA_gff  mitos.result.update.tRNA.gff  -mfannot_tRNA_gff  mfannot.update.tRNA.gff  -tRNAscan_out  tRNAscan-SE.out.update.txt  -output_prefix  tRNA 
version: 1.0
Options:
        -mitos_tRNA_gff <file>                 mitos predicted tRNA, such as 'mitos.result.update.tRNA.gff'
        -mfannot_tRNA_gff <file>               mfanno predicted tRNA, such as 'mfannot.update.tRNA.gff'
        -tRNAscan_out <file>                   tRNAscan predicted tRNA, such as 'tRNAscan-SE.out.update.txt'
        -output_prefix <strings>               output prefix, such as 'tRNA'
        -help                                  print help information

USAGE
############################## usage end #######################################

if ($help || !(defined $mitos_tRNA_gff) || !(defined $mfannot_tRNA_gff) || !(defined $tRNAscan_out) || !(defined $output_prefix))
{
	print $usage;
	exit;
}

my @months = qw(Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec);
my ($timesecond, $timeminute, $timehour, $timedaymonths, $timemonth, $timeyear, $timedayweek, $timedayYear, $timedaylightsavings) = localtime();
$timeyear+=1900;

print "[$timehour\:$timeminute\:$timesecond\, $months[$timemonth] $timedaymonths\, $timeyear] ";
print "Start to run the \'$0\' program ...\n";

# mitos.result.update.tRNA.gff
# CP084950.1	mitfi	tRNA	1672	1740	.	+	.	ID=tRNA:trnF_1;Name=trnF_1

# mfannot.update.gene.gff
# CP084950.1	MFAnnot	tRNA	1672	1740	.	+	.	ID=tRNA:trnF(aaa);Name=trnF(aaa)

# tRNAscan-SE.out.update.txt
# CP084950.1	1	1672 	1740 	Phe	AAA	0	0	28.1	

################################################################################
# Reading 'tRNAscan_out'
################################################################################

print "Reading '$tRNAscan_out' file ...\n";

my %pos;
my %tRNA;
open SE,"$tRNAscan_out" || die "Cannot open the file '$tRNAscan_out'.\n";
while(<SE>)
{
	chomp;
	my @sp=split(/\t/,$_);
	if($_ ne "" && !($_=~/^Sequence/ || $_=~/^Name/ || $_=~/^--------/))
	{
		my $strand="+";
		if($sp[3]<$sp[2])
		{
			$strand="-";
		}
		if($strand eq "+")
		{
			for(my $i=$sp[2];$i<=$sp[3];$i++)
			{
				$pos{$i}=1;
			}
			$tRNA{$sp[0]}{$sp[2]}="$sp[0]\ttRNAscan-SE\ttRNA\t$sp[2]\t$sp[3]\t\.\t+\t\.\tID=tRNA:$sp[4]\($sp[5]\);Name=$sp[4]\($sp[5]\)";
		}
		else
		{
			for(my $i=$sp[3];$i<=$sp[2];$i++)
			{
				$pos{$i}=1;
			}
			$tRNA{$sp[0]}{$sp[2]}="$sp[0]\ttRNAscan-SE\ttRNA\t$sp[3]\t$sp[2]\t\.\t-\t\.\tID=tRNA:$sp[4]\($sp[5]\);Name=$sp[4]\($sp[5]\)";
		}
	}
}
close SE;

################################################################################
# Reading 'mitos_tRNA_gff'
################################################################################

print "Reading '$mitos_tRNA_gff' file ...\n";

open MI,"$mitos_tRNA_gff" || die "Cannot open the file '$mitos_tRNA_gff'.\n";
while(<MI>)
{
	chomp;
	my @sp=split(/\t/,$_);
	my $tag=0;
	if($_ ne "")
	{
		for(my $i=$sp[3];$i<=$sp[4];$i++)
		{
			if(defined $pos{$i})
			{
				$tag=1;
				last;
			}
		}
		if($tag==0)
		{
			for(my $i=$sp[3];$i<=$sp[4];$i++)
			{
				$pos{$i}=1;
			}
			$tRNA{$sp[0]}{$sp[3]}=$_;
		}
	}
}
close MI;

################################################################################
# Reading 'mfannot_tRNA_gff'
################################################################################

print "Reading '$mfannot_tRNA_gff' file ...\n";

open MF,"$mfannot_tRNA_gff" || die "Cannot open the file '$mfannot_tRNA_gff'.\n";
while(<MF>)
{
	chomp;
	my @sp=split(/\t/,$_);
	my $tag=0;
	if($_ ne "")
	{
		for(my $i=$sp[3];$i<=$sp[4];$i++)
		{
			if(defined $pos{$i})
			{
				$tag=1;
				last;
			}
		}
		if($tag==0)
		{
			for(my $i=$sp[3];$i<=$sp[4];$i++)
			{
				$pos{$i}=1;
			}
			$tRNA{$sp[0]}{$sp[3]}=$_;
		}
	}
}
close MF;

open OA,">$output_prefix\_ori.gff";
foreach my $sid(sort keys %tRNA)
{
	foreach my $pos(sort {$a<=>$b} keys %{$tRNA{$sid}})
	{
		print OA "$tRNA{$sid}{$pos}\n";
	}
}
close OA;

# CP084950.1	tRNAscan-SE	tRNA	1672 	1740 	.	+	.	ID=tRNA:Phe(AAA);Name=Phe(AAA)
my %num;
open OB,">$output_prefix\_result.gff";
open G,"$output_prefix\_ori.gff" || die "Cannot open the file '$output_prefix\_ori.gff'.\n";
while(<G>)
{
	chomp;
	my @sp=split(/\t/,$_);
	if($_ ne "")
	{
		my $id=(split(/\;/,$sp[8]))[0];
		$id=~s/ID=tRNA://;
		if(!(defined $num{$id}))
		{
			print OB "$sp[0]\t$sp[1]\t$sp[2]\t$sp[3]\t$sp[4]\t$sp[5]\t$sp[6]\t$sp[7]\tID=tRNA\:$id\_1\;Name=$id\_1\n";
			$num{$id}++;
		}
		else
		{
			my $fre=$num{$id};
			$fre++;
			print OB "$sp[0]\t$sp[1]\t$sp[2]\t$sp[3]\t$sp[4]\t$sp[5]\t$sp[6]\t$sp[7]\tID=tRNA\:$id\_$fre\;Name=$id\_$fre\n";
			$num{$id}++;
		}
	}
}
close G;
close OB;

($timesecond, $timeminute, $timehour, $timedaymonths, $timemonth, $timeyear, $timedayweek, $timedayYear, $timedaylightsavings) = localtime();
$timeyear+=1900;
print "[$timehour\:$timeminute\:$timesecond\, $months[$timemonth] $timedaymonths\, $timeyear] ";
print "End of running the \'$0\' program.\n";

__END__
