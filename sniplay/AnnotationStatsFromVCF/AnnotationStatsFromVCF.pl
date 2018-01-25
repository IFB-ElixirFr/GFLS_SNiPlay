#!/usr/bin/perl

use strict;
use Switch;
use Getopt::Long;
use Bio::SeqIO;

my $usage = qq~Usage:$0 <args> [<opts>]
where <args> are:
    -v, --vcf           <VCF input>
    -o, --out           <output>
    -s, --step          <step (in bp)>
~;
$usage .= "\n";

my ($vcf,$out,$step,$min_depth);
$min_depth = 5;

GetOptions(
	"vcf=s"      => \$vcf,
	"out=s"      => \$out,
	"min_depth=s"=> \$min_depth,
	"step=s"     => \$step,
);


die $usage
  if ( !$vcf || !$step || !$out );

if ($step =~/^(\d+)\s*$/){
	$step = $1;
}
else{
	die "Error: step size must be an integer\n";
}


my $VCFTOOLS_EXE = "vcftools"; 
 
my %max;
my %counts_ns;
my %counts_s;
my $nb_gene = 0;
my $nb_intergenic = 0;
my $nb_exon = 0;
my $nb_intron = 0;
my $nb_UTR = 0;
my $nb_syn = 0;
my $nb_nsyn = 0;
if ($vcf =~/\.bcf/){
	system("$VCFTOOLS_EXE --bcf $vcf --recode --recode-INFO-all --out $out");
	$vcf = "$out.recode.vcf";
}
open(my $VCF,$vcf);
while(<$VCF>)
{
	my @infos = split(/\t/,$_);
	if (scalar @infos > 8 && !/#CHROM/)
	{
		my $chrom = $infos[0];
		my $position = $infos[1];
		my $id = $infos[2];
		
		
		my $classe_position = int($position/$step);
		if (/=NON_SYNONYMOUS_CODING/){
			$counts_ns{$chrom}{$classe_position}++;
			$nb_nsyn++;
		}
		if (/=SYNONYMOUS_CODING/){
			$counts_s{$chrom}{$classe_position}++;
			$nb_syn++;
		}
		if (/INTERGENIC/){
			$nb_intergenic++;
		}
		elsif (/EFF=/){
			$nb_gene++;
		}
		if (/CODING/){
			}
		elsif (/UTR/){
			$nb_UTR++;
		}
		elsif (/INTRON/){
			$nb_intron++;
		}
		$max{$chrom} = $classe_position;
	}
}
my $nb_exon = $nb_gene - $nb_intron - $nb_UTR;

open(my $OUT,">$out");
#print $OUT "Chrom	Bin	Nb synonymous SNPs	Nb non-synonymous SNPs	dN/dS ratio\n";
print $OUT "Chrom	Bin	dN/dS ratio\n";
foreach my $chrom(sort keys(%counts_s))
{
	my $maximum = $max{$chrom};
	for(my $i=1;$i<=$maximum;$i++)
	{
		my $classe_position = $i;
		my $nb_s = 0;
		my $nb_ns = 0;
		my $ratio = 0;
		if ($counts_s{$chrom}{$classe_position}){$nb_s = $counts_s{$chrom}{$classe_position};}
		if ($counts_ns{$chrom}{$classe_position}){$nb_ns = $counts_ns{$chrom}{$classe_position};}
		if ($nb_s){$ratio = $nb_ns/$nb_s;}
		my $bin = $classe_position * $step;
		#print $OUT "$chrom	$classe_position	$nb_s	$nb_ns	$ratio\n";
		print $OUT "$chrom	$bin	$ratio\n";
	}
}
close($OUT);



open(my $A2, ">$out.location");
print $A2 "Intergenic	$nb_intergenic	Intergenic:$nb_intergenic\n";
print $A2 "Genic	$nb_gene	Exon:$nb_exon	Intron:$nb_intron	UTR:$nb_UTR\n";
close($A2);

open(my $A2, ">$out.effect");
print $A2 "Intron	$nb_intron	Intron:$nb_intron\n";
print $A2 "UTR	$nb_UTR	UTR:$nb_UTR\n";
print $A2 "Exon	$nb_exon	Synonym:$nb_syn	Non-syn:$nb_nsyn\n";
close($A2);

