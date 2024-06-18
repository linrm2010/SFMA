#! /usr/bin/perl -w
################################################################################
# mitogenome_gff_product_ncRNA2tbl.pl version 1.0
# Copyright (C) Runmao Lin, Tong Liu, Xiaoting Wang, Fanxing Yang, Zhiyin Wang, 2023
# Contact (E-mail): linrunmao@hainanu.cn
#
# This program is provided under the terms of a personal license to the recipient and
# may only be used for the recipient's own research at an academic insititution.
#
# For using this program in a company or for commercial purposes, a commercial license
# is required.
#
# The purpose of this program is to generate tbl file for uploading genbank to NCBI.
################################################################################

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;

my $gene_gff;
my $gene_product;
my $tRNA_gff;
my $rRNA_gff;
my $dbname;
my $output_file;
my $help;

GetOptions
(
	"gene_gff=s" => \$gene_gff,                             # string
	"gene_product=s" => \$gene_product,                     # string
	"tRNA_gff=s" => \$tRNA_gff,                             # string
	"rRNA_gff=s" => \$rRNA_gff,                             # string
	"dbname=s" => \$dbname,                                 # string
	"output_file=s" => \$output_file,                       # string
	"help" => \$help                                        # flag
);

################################################################################
# usage
############################# usage begin ######################################
my $usage= << "USAGE";

Example: perl $0  -gene_gff  gene.gff  -gene_product  gene_product.txt  -tRNA_gff  tRNA.gff  -rRNA_gff rRNA.gff  -dbname HNU  -output_file  genome.tble 
version: 1.0
Options:
        -gene_gff <file>                 gene gff file, such as 'gene.gff'
        -gene_product <file>             gene product, such as 'gene_product.txt'
        -tRNA_gff <file>                 tRNA gff, such as 'tRNA.gff'
        -rRNA_gff <file>                 rRNA gff, such as 'rRNA.gff'
        -dbname <strings>                such as 'HNU'
        -output_file <file>              output file, such as 'genome.tbl'
        -help                            print help information

USAGE
############################## usage end #######################################

if ($help || !(defined $gene_gff) || !(defined $gene_product) || !(defined $tRNA_gff) || !(defined $rRNA_gff) || !(defined $dbname) || !(defined $output_file))
{
	print $usage;
	exit;
}

my @months = qw(Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec);
my ($timesecond, $timeminute, $timehour, $timedaymonths, $timemonth, $timeyear, $timedayweek, $timedayYear, $timedaylightsavings) = localtime();
$timeyear+=1900;

print "[$timehour\:$timeminute\:$timesecond\, $months[$timemonth] $timedaymonths\, $timeyear] ";
print "Start to run the \'$0\' program ...\n";

### incomplete gene ###
# >Feature contig
# <pos1	>pos2	gene
# 			locus_tag	geneid
# <posa	posb	mRNA
# ...
# posa	>posb	
# 			product	...
# 			note	UniProt
# 			protein_id	gnl|SCAUAG1IA|geneid
# 			transcript_id	gnl|SCAUAG1IA|mrna.geneid
# posa	posb	CDS
# ...
# posa	posb	
# 			product	...
# 			note	UniProt
# 			protein_id	gnl|SCAUAG1IA|geneid
# 			transcript_id	gnl|SCAUAG1IA|mrna.geneid
# 			go_process	...
# 			go_component	...
# 			go_function	...
#
### complete gene ###
# >Feature contig
# pos1	pos2	gene
# 			locus_tag	geneid
# posa	posb	mRNA
# ...
# posa	posb	
# 			product	...
# 			note	UniProt
# 			protein_id	gnl|SCAUAG1IA|geneid
# 			transcript_id	gnl|SCAUAG1IA|mrna.geneid
# posa	posb	CDS
# ...
# posa	posb	
# 			product	...
# 			note	UniProt
# 			protein_id	gnl|SCAUAG1IA|geneid
# 			transcript_id	gnl|SCAUAG1IA|mrna.geneid
# 			go_process	...
# 			go_component	...
# 			go_function	...
#
### ncRNA ###
# pos	pos	gene
# 			locus_tag	geneid
# 			pseudo
# pos	pos	tRNA/rRNA
# 			product	Leu/18S/...
# 

my %contig;

# mRNA product file
print "Reading product ...\n";
my %product;
my %note;
open P,"$gene_product" || die "Cannot open the file '$gene_product'.\n";
while(<P>)
{
	chomp;
	my @sp=split(/\t/,$_);
	if($_ ne "" && !($_=~/^GeneID/))
	{
		$product{$sp[0]}=$sp[1];
		if(defined $sp[2] && !($sp[2]=~/^PF/ || $sp[2] eq ""))
#		if(defined $sp[2] && $sp[2] ne "")
		{
			$note{$sp[0]}=$sp[3];
		}
		else
		{
			$note{$sp[0]}="";
		}
	}
}
close P;

# mRNA gff file
print "Reading GFF ...\n";
my $geneid="";
my %mRNA; # "AG1IA_00001"=>"
			# <pos1	>pos2	gene
			# 			locus_tag	geneid
			# <posa	posb	mRNA
			# ...
			# posa	>posb	
			# 			product	...
			# 			note	UniProt
			# 			protein_id	gnl|SCAUAG1IA|geneid
			# 			transcript_id	gnl|SCAUAG1IA|mrna.geneid
			# posa	posb	CDS
			# ...
			# posa	posb	
			# 			product	...
			# 			note	UniProt
			# 			protein_id	gnl|SCAUAG1IA|geneid
			# 			transcript_id	gnl|SCAUAG1IA|mrna.geneid
my %genebegin; # $genebegin{"contig1"}{pos1}="AG1IA_00001"
open CG,"$gene_gff" || die "Cannot open the file '$gene_gff'.\n";
while(<CG>)
{
	chomp;
	my @sp=split(/\t/,$_);
	if($_=~/\tmRNA\t/)
	{
		$geneid=(split(/\;/,$sp[8]))[0];
		$geneid=~s/ID=mRNA\://;
		my $tcid=$sp[0];
		$contig{$tcid}=1;
		my $contigbegin=0;
		my $contigend=0;
#		print "***$geneid***\n";
#		my $utr3=0;
#		my $utr5=0;
		my @pos;
		while(<CG>)
		{
			chomp;
			my @gsp=split(/\t/,$_);
			if($_=~/\tCDS\t/)
			{
				$contigbegin=$gsp[3] if($contigbegin==0);
				$contigend=$gsp[4];
				if($sp[6] eq "+")
				{
					push(@pos,"$gsp[3]\t$gsp[4]");
				}
				else
				{
					push(@pos,"$gsp[4]\t$gsp[3]");
				}
			}
#			$utr3++ if($_=~/three_prime_UTR/);
#			$utr5++ if($_=~/five_prime_UTR/);
#			last if($utr3>0 && $utr5>0);
			last if($_=~/\tgene\t/);
		}
		if($sp[6] eq "+")
		{
			$mRNA{$geneid}="\<$contigbegin\t\>$contigend\tgene\t\t\n";
			$mRNA{$geneid}.="\t\t\tlocus_tag\t$geneid\n";
			if(scalar @pos==1)
			{
				$mRNA{$geneid}.="\<$contigbegin\t\>$contigend\tmRNA\t\t\n";
				$mRNA{$geneid}.="\t\t\tproduct\t$product{$geneid}\n";
				$mRNA{$geneid}.="\t\t\tnote\t$note{$geneid}\n" if(defined $note{$geneid} && $note{$geneid} ne "");
				$mRNA{$geneid}.="\t\t\tprotein_id\tgnl\|$dbname\|$geneid\n"; # protein_id	gnl|SCAUAG1IA|geneid
				$mRNA{$geneid}.="\t\t\ttranscript_id\tgnl\|$dbname\|mrna\.$geneid\n";   #transcript_id	gnl|SCAUAG1IA|mrna.geneid
				$mRNA{$geneid}.="$contigbegin\t$contigend\tCDS\t\t\n";
				$mRNA{$geneid}.="\t\t\tproduct\t$product{$geneid}\n";
				$mRNA{$geneid}.="\t\t\tnote\t$note{$geneid}\n" if(defined $note{$geneid} && $note{$geneid} ne "");
				$mRNA{$geneid}.="\t\t\tprotein_id\tgnl\|$dbname\|$geneid\n"; # protein_id	gnl|SCAUAG1IA|geneid
				$mRNA{$geneid}.="\t\t\ttranscript_id\tgnl\|$dbname\|mrna\.$geneid\n";   #transcript_id	gnl|SCAUAG1IA|mrna.geneid
			}
			else
			{
				$mRNA{$geneid}.="\<$pos[0]\tmRNA\t\t\n";
				for(my $i=1;$i<@pos-1;$i++)
				{
					$mRNA{$geneid}.="$pos[$i]\t\t\t\n";
				}
				my @spos=split(/\t/,$pos[-1]);
				$mRNA{$geneid}.="$spos[0]\t\>$spos[1]\t\t\t\n";
				$mRNA{$geneid}.="\t\t\tproduct\t$product{$geneid}\n";
				$mRNA{$geneid}.="\t\t\tnote\t$note{$geneid}\n" if(defined $note{$geneid} && $note{$geneid} ne "");
				$mRNA{$geneid}.="\t\t\tprotein_id\tgnl\|$dbname\|$geneid\n"; # protein_id	gnl|SCAUAG1IA|geneid
				$mRNA{$geneid}.="\t\t\ttranscript_id\tgnl\|$dbname\|mrna\.$geneid\n";   #transcript_id	gnl|SCAUAG1IA|mrna.geneid
				$mRNA{$geneid}.="$pos[0]\tCDS\t\t\n";
				for(my $i=1;$i<@pos;$i++)
				{
					$mRNA{$geneid}.="$pos[$i]\t\t\t\n";
				}
				$mRNA{$geneid}.="\t\t\tproduct\t$product{$geneid}\n";
				$mRNA{$geneid}.="\t\t\tnote\t$note{$geneid}\n" if(defined $note{$geneid} && $note{$geneid} ne "");
				$mRNA{$geneid}.="\t\t\tprotein_id\tgnl\|$dbname\|$geneid\n"; # protein_id	gnl|SCAUAG1IA|geneid
				$mRNA{$geneid}.="\t\t\ttranscript_id\tgnl\|$dbname\|mrna\.$geneid\n";   #transcript_id	gnl|SCAUAG1IA|mrna.geneid
			}
		}
		else
		{
			$mRNA{$geneid}="\<$contigend\t\>$contigbegin\tgene\t\t\n";
			$mRNA{$geneid}.="\t\t\tlocus_tag\t$geneid\n";
			if(scalar @pos==1)
			{
				$mRNA{$geneid}.="\<$contigend\t\>$contigbegin\tmRNA\t\t\n";
				$mRNA{$geneid}.="\t\t\tproduct\t$product{$geneid}\n";
				$mRNA{$geneid}.="\t\t\tnote\t$note{$geneid}\n" if(defined $note{$geneid} && $note{$geneid} ne "");
				$mRNA{$geneid}.="\t\t\tprotein_id\tgnl\|$dbname\|$geneid\n"; # protein_id	gnl|SCAUAG1IA|geneid
				$mRNA{$geneid}.="\t\t\ttranscript_id\tgnl\|$dbname\|mrna\.$geneid\n";   #transcript_id	gnl|SCAUAG1IA|mrna.geneid
				$mRNA{$geneid}.="$contigend\t$contigbegin\tCDS\t\t\n";
				$mRNA{$geneid}.="\t\t\tproduct\t$product{$geneid}\n";
				$mRNA{$geneid}.="\t\t\tnote\t$note{$geneid}\n" if(defined $note{$geneid} && $note{$geneid} ne "");
				$mRNA{$geneid}.="\t\t\tprotein_id\tgnl\|$dbname\|$geneid\n"; # protein_id	gnl|SCAUAG1IA|geneid
				$mRNA{$geneid}.="\t\t\ttranscript_id\tgnl\|$dbname\|mrna\.$geneid\n";   #transcript_id	gnl|SCAUAG1IA|mrna.geneid
			}
			else
			{
				$mRNA{$geneid}.="\<$pos[-1]\tmRNA\t\t\n";
				for(my $i=@pos-2;$i>0;$i--)
				{
					$mRNA{$geneid}.="$pos[$i]\t\t\t\n";
				}
				my @spos=split(/\t/,$pos[0]);
				$mRNA{$geneid}.="$spos[0]\t\>$spos[1]\t\t\t\n";
				$mRNA{$geneid}.="\t\t\tproduct\t$product{$geneid}\n";
				$mRNA{$geneid}.="\t\t\tnote\t$note{$geneid}\n" if(defined $note{$geneid} && $note{$geneid} ne "");
				$mRNA{$geneid}.="\t\t\tprotein_id\tgnl\|$dbname\|$geneid\n"; # protein_id	gnl|SCAUAG1IA|geneid
				$mRNA{$geneid}.="\t\t\ttranscript_id\tgnl\|$dbname\|mrna\.$geneid\n";   #transcript_id	gnl|SCAUAG1IA|mrna.geneid
				$mRNA{$geneid}.="$pos[-1]\tCDS\t\t\n";
				for(my $i=@pos-2;$i>=0;$i--)
				{
					$mRNA{$geneid}.="$pos[$i]\t\t\t\n";
				}
				$mRNA{$geneid}.="\t\t\tproduct\t$product{$geneid}\n";
				$mRNA{$geneid}.="\t\t\tnote\t$note{$geneid}\n" if(defined $note{$geneid} && $note{$geneid} ne "");
				$mRNA{$geneid}.="\t\t\tprotein_id\tgnl\|$dbname\|$geneid\n"; # protein_id	gnl|SCAUAG1IA|geneid
				$mRNA{$geneid}.="\t\t\ttranscript_id\tgnl\|$dbname\|mrna\.$geneid\n";   #transcript_id	gnl|SCAUAG1IA|mrna.geneid
			}
		}
		if($contigbegin<$contigend)
		{
			$genebegin{$sp[0]}{$contigbegin}=$geneid;
		}
		else
		{
			$genebegin{$sp[0]}{$contigend}=$geneid;
		}
	}
}
close CG;

# ncRNA
my %ncRNA; # "AG1IA_10490"=>
			# pos	pos	gene
			# 			locus_tag	geneid
			# 			pseudo
			# pos	pos	tRNA/rRNA
			# 			product	Leu/18S/...
# tRNA
print "Reading tRNA ...\n";
open TR,"$tRNA_gff" || die "Cannot open the file '$tRNA_gff'.\n";
while(<TR>)
{
	chomp;
	my @sp=split(/\t/,$_);
	if($_ ne "")
	{
		my $tcid=$sp[0];
		my $tRNA_id=(split(/\;/,$sp[8]))[0];
		$tRNA_id=~s/ID=tRNA://;
		$contig{$tcid}=1;
		$sp[3]=~s/\s+//g;
		$sp[4]=~s/\s+//g;
		$ncRNA{$tRNA_id}="$sp[3]\t$sp[4]\tgene\t\t\n";
		$ncRNA{$tRNA_id}.="\t\t\tlocus_tag\t$tRNA_id\n";
#		$ncRNA{$tRNA_id}.="\t\t\tpseudo\t\n" if($sp[8] eq "Y" || $sp[8] eq "y");
		$ncRNA{$tRNA_id}.="$sp[3]\t$sp[4]\ttRNA\t\t\n";
		$ncRNA{$tRNA_id}.="\t\t\tproduct\t$tRNA_id\n";
		$genebegin{$sp[0]}{$sp[3]}=$tRNA_id;
	}
}
close TR;
# rRNA
my %rRNA_infor;
my %rRNA_infor_intron;
print "Reading rRNA ...\n";
open RR,"$rRNA_gff" || die "Cannot open the file '$rRNA_gff'.\n";
while(<RR>)
{
	chomp;
	my @sp=split(/\t/,$_);
	if($_=~/\trRNA\t/)
	{
		my $tcid=$sp[0];
		$contig{$tcid}=1;
		my $rRNA_id=(split(/\;/,$sp[8]))[0];
		$rRNA_id=~s/ID=rRNA://;
#		$ncRNA{$rRNA_id}="$sp[3]\t$sp[4]\tgene\t\t\n";
#		$ncRNA{$rRNA_id}.="\t\t\tlocus_tag\t$rRNA_id\n";
#		$ncRNA{$rRNA_id}.="$sp[3]\t$sp[4]\trRNA\t\t\n";
#		$ncRNA{$rRNA_id}.="\t\t\tproduct\t$rRNA_id\n";
		$rRNA_infor{$rRNA_id}="$_";
		$genebegin{$sp[0]}{$sp[3]}=$rRNA_id;
	}
	elsif($_=~/\tintron\t/)
	{
		my $rRNA_id=(split(/\;/,$sp[8]))[0];
		$rRNA_id=~s/ID=rRNA://;
		$rRNA_infor_intron{$rRNA_id}.="$_\n";
	}
}
close RR;

foreach my $rRNA_id(keys %rRNA_infor)
{
	my @rsp=split(/\t/,$rRNA_infor{$rRNA_id});
	if(!(defined $rRNA_infor_intron{$rRNA_id}))
	{
		$ncRNA{$rRNA_id}="$rsp[3]\t$rsp[4]\tgene\t\t\n";
		$ncRNA{$rRNA_id}.="\t\t\tlocus_tag\t$rRNA_id\n";
		$ncRNA{$rRNA_id}.="$rsp[3]\t$rsp[4]\trRNA\t\t\n";
		$ncRNA{$rRNA_id}.="\t\t\tproduct\t$rRNA_id\n";
	}
	else
	{
		my %intron_list;
		my @isp=split(/\n/,$rRNA_infor_intron{$rRNA_id});
		for(my $i=0;$i<@isp;$i++)
		{
			my @sisp=split(/\t/,$isp[$i]);
			$intron_list{$sisp[3]}="$sisp[3] $sisp[4]";
		}
		my $rRNA_begin=$rsp[3];
		my $rRNA_end=$rsp[4];
		my $tag_first=0;
		$ncRNA{$rRNA_id}="$rsp[3]\t$rsp[4]\tgene\t\t\n";
		$ncRNA{$rRNA_id}.="\t\t\tlocus_tag\t$rRNA_id\n";
		foreach my $ipos(sort{$a<=>$b} keys %intron_list)
		{
			my @sppos=split(/\s+/,$intron_list{$ipos});
			if($tag_first==0)
			{
				$ncRNA{$rRNA_id}.="$rRNA_begin\t".($sppos[0]-1)."\trRNA\t\t\n";
				$tag_first=$sppos[1]+1;
				$rRNA_begin=$sppos[1]+1;
			}
			else
			{
				$ncRNA{$rRNA_id}.="$rRNA_begin\t".($sppos[0]-1)."\t\t\t\n";
				$tag_first=$sppos[1]+1;
				$rRNA_begin=$sppos[1]+1;
			}
		}
		$ncRNA{$rRNA_id}.="$rRNA_begin\t$rRNA_end\t\t\t\n";
		$ncRNA{$rRNA_id}.="\t\t\tproduct\t$rRNA_id\n";
	}
}

# output
print "Preparing for output ...\n";
open OUT,">$output_file";
foreach my $idnum(sort keys %contig)
{
	my $id="$idnum";
	print OUT "\>Features $id\n";
	foreach my $loc(sort{$a<=>$b} keys %{$genebegin{$id}})
	{
		print OUT "$mRNA{$genebegin{$id}{$loc}}" if(defined $mRNA{$genebegin{$id}{$loc}});
		print OUT "$ncRNA{$genebegin{$id}{$loc}}" if(defined $ncRNA{$genebegin{$id}{$loc}});
	}
}
close OUT;

($timesecond, $timeminute, $timehour, $timedaymonths, $timemonth, $timeyear, $timedayweek, $timedayYear, $timedaylightsavings) = localtime();
$timeyear+=1900;
print "[$timehour\:$timeminute\:$timesecond\, $months[$timemonth] $timedaymonths\, $timeyear] ";
print "End of running the \'$0\' program.\n";

__END__
