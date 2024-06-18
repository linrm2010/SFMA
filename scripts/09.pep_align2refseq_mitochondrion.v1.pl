#! /usr/bin/perl -w
################################################################################
# pep_align2refseq_mitochondrion.pl version 1.0
# Copyright (C) Runmao Lin, Tong Liu, Xiaoting Wang, Fanxing Yang, Zhiyin Wang, 2023
# Contact (E-mail): linrunmao@hainanu.cn
#
# This program is provided under the terms of a personal license to the recipient and
# may only be used for the recipient's own research at an academic insititution.
#
# For using this program in a company or for commercial purposes, a commercial license
# is required.
#
# The purpose of this program is to display alignment between predicted genes and refseq genes.
################################################################################

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;

my $gene_pep;
my $refseq_mitochondrion_pep;
my $blast_tab;
my $output_prefix;
my $help;

GetOptions
(
	"gene_pep=s" => \$gene_pep,                                         # string
	"refseq_mitochondrion_pep=s" => \$refseq_mitochondrion_pep,         # string
	"blast_tab=s" => \$blast_tab,                                       # string
	"output_prefix=s" => \$output_prefix,                               # string
	"help" => \$help                                                    # flag
);

################################################################################
# usage
############################# usage begin ######################################
my $usage= << "USAGE";

Example: perl $0 -gene_pep  gene.pep  -refseq_mitochondrion_pep  refseq_mitochondrion.pep  -blast_tab  blast.tab  -output_prefix  mitos 
version: 1.0
Options:
        -gene_pep <file>                       gene pep file, such as 'gene.pep'
        -refseq_mitochondrion_pep <file>       mitochondrial genes in refseq database, such as 'refseq_mitochondrion.pep'
        -blast_tab <file>                      alignment tab file by BLAST, such as 'blast.tab'
        -output_prefix <strings>               such as 'mitos'
        -help                                  print help information

USAGE
############################## usage end #######################################

if ($help || !(defined $gene_pep) || !(defined $refseq_mitochondrion_pep) || !(defined $blast_tab) || !(defined $output_prefix))
{
	print $usage;
	exit;
}

my @months = qw(Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec);
my ($timesecond, $timeminute, $timehour, $timedaymonths, $timemonth, $timeyear, $timedayweek, $timedayYear, $timedaylightsavings) = localtime();
$timeyear+=1900;

print "[$timehour\:$timeminute\:$timesecond\, $months[$timemonth] $timedaymonths\, $timeyear] ";
print "Start to run the \'$0\' program ...\n";

my $gene_len=fasta_seq_length($gene_pep);
my $refseq_len=fasta_seq_length($refseq_mitochondrion_pep);

# Query id.1	Subject id.2	identity.3	alignment length.4	mismatches.5	gap openings.6	q. start.7	q. end.8	s. start.9	s. end.10	e-value.11	bit score.12
# atp6	YP_009254348.1	98.46	260	4	0	1	260	1	260	0.0	 499

my %select_num; # used top 20 alignment for drawing maps
my %select_infor;
open F,"$blast_tab" || die "Cannot open the file '$blast_tab'.\n";
while(<F>)
{
	chomp;
	my @sp=split(/\t/,$_);
	if($_ ne "" && (!(defined $select_num{$sp[0]}) || (defined $select_num{$sp[0]} && $select_num{$sp[0]}<20)))
	{
		$select_num{$sp[0]}++;
		$select_infor{$sp[0]}.="$_\n";
	}
}
close F;

foreach my $id(sort keys %select_infor)
{
	open OA,">$output_prefix\_$id.svg";
	print OA "<?xml  version=\"1.0\"  standalone=\"no\"  ?>\n";
	print OA "<!DOCTYPE  svg  PUBLIC  \"-//W3C//DTD  SVG  1.0//EN\" \"http://www.w3.org/TR/2001/REC-SVG-20010904/DTD/svg10.dtd\">\n";
	print OA "<svg  width=\"655\"  height=\"378\"  xmlns=\"http://www.w3.org/2000/svg\">\n";
	print OA "<text  x=\"260\"  y=\"20\"  font-family=\"Arial\" font-size=\"8.75\">BLAST alignment of $id against other genes</text>\n";
	print OA "<line x1=\"450\" y1=\"30\" x2=\"550\" y2=\"30\" stroke=\"red\"  stroke-width=\"1\" />\n";
	print OA "<text  x=\"550\"  y=\"30\"  font-family=\"Arial\" font-size=\"8.75\">matched regions</text>\n";
	my @sp=split(/\n/,$select_infor{$id});
	my @asp=split(/\t/,$sp[0]);
	print OA "<rect x=\"10\"  y=\"50\"  width=\"250\"  height=\"10\"  fill=\"blue\" stroke=\"none\"/>\n";
	print OA "<text  x=\"100\"  y=\"45\"  font-family=\"Arial\" font-size=\"8.75\">$id</text>\n";
	print OA "<text  x=\"10\"  y=\"50\"  font-family=\"Arial\" font-size=\"8.75\">1</text>\n";
	print OA "<text  x=\"260\"  y=\"50\"  font-family=\"Arial\" font-size=\"8.75\">$gene_len->{$id}</text>\n";
	print OA "<text  x=\"400\"  y=\"45\"  font-family=\"Arial\" font-size=\"8.75\">Subject</text>\n";
	print OA "<rect x=\"300\"  y=\"50\"  width=\"250\"  height=\"10\"  fill=\"blue\" stroke=\"none\"/>\n";
	for(my $i=0;$i<@sp;$i++)
	{
		my $query_x_begin=10;
		my $query_y_begin=75+$i*15;
		my $subject_x_begin=300;
		my $subject_y_begin=$query_y_begin;
		my @psp=split(/\t/,$sp[$i]);
		my $query_each_base=250/$gene_len->{$id};
		my $subject_each_base=250/$refseq_len->{$psp[1]};
		print OA "<line x1=\"10\" y1=\"$query_y_begin\" x2=\"260\" y2=\"$query_y_begin\" stroke=\"blue\"  stroke-width=\"1\" />\n";
		print OA "<line x1=\"".($query_x_begin+($psp[6]-1)*$query_each_base)."\" y1=\"".($query_y_begin+3)."\" x2=\"".($query_x_begin+($psp[7]-1)*$query_each_base)."\" y2=\"".($query_y_begin+3)."\" stroke=\"red\"  stroke-width=\"1\" />\n";
		print OA "<line x1=\"300\" y1=\"$subject_y_begin\" x2=\"550\" y2=\"$subject_y_begin\" stroke=\"blue\"  stroke-width=\"1\" />\n";
		print OA "<text  x=\"300\"  y=\"$subject_y_begin\"  font-family=\"Arial\" font-size=\"8.75\">1</text>\n";
		print OA "<text  x=\"550\"  y=\"$subject_y_begin\"  font-family=\"Arial\" font-size=\"8.75\">$refseq_len->{$psp[1]}</text>\n";
		print OA "<text  x=\"580\"  y=\"$subject_y_begin\"  font-family=\"Arial\" font-size=\"8.75\">$psp[1]</text>\n";
		print OA "<line x1=\"".($subject_x_begin+($psp[8]-1)*$subject_each_base)."\" y1=\"".($subject_y_begin+3)."\" x2=\"".($subject_x_begin+($psp[9]-1)*$subject_each_base)."\" y2=\"".($subject_y_begin+3)."\" stroke=\"red\"  stroke-width=\"1\" />\n";
	}
	print OA "<\/svg>\n";
	close OA;
#	system(" convert -density 300 $ARGV[3] $ARGV[3]\.png ");
}


sub fasta_seq_length
{
	my $file=shift;
	my $len;
	my $seq;
	my $seq_id="";
	open S,"$file" || die "Cannot open the file '$file'.\n";
	while(<S>)
	{
		chomp;
		if($_=~/>/)
		{
			$seq_id=(split(/\s+/,$_))[0];
			$seq_id=~s/>//;
		}
		elsif($_ ne "")
		{
			$seq->{$seq_id}.="$_";
			$seq->{$seq_id}=~s/\s+//;
		}
	}
	close S;
	foreach my $id(keys %{$seq})
	{
		$len->{$id}=length $seq->{$id};
	}
	return $len;
}

($timesecond, $timeminute, $timehour, $timedaymonths, $timemonth, $timeyear, $timedayweek, $timedayYear, $timedaylightsavings) = localtime();
$timeyear+=1900;
print "[$timehour\:$timeminute\:$timesecond\, $months[$timemonth] $timedaymonths\, $timeyear] ";
print "End of running the \'$0\' program.\n";

__END__
