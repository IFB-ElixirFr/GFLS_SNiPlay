#!/usr/bin/perl

use strict;
use Getopt::Long;
use Bio::SeqIO;

my $usage = qq~Usage:$0 <args> [<opts>]
where <args> are:
    -i, --input         <VCF/Hapmap input>
    -o, --out           <output in tabular format>
    -s, --step          <step (in bp)>
~;
$usage .= "\n";

my ($input,$out,$step);

GetOptions(
	"input=s"    => \$input,
	"out=s"      => \$out,
	"step=s"     => \$step,
);


die $usage
  if ( !$input || !$step || !$out );
  
my $max_chr_num = 100;

my %counts;
my %counts_by_ind;

my $VCFTOOLS_EXE = "vcftools";

my $is_vcf = `head -4000 $input | grep -c CHROM`;
my $is_bcf = 0;
if ($input =~/\.bcf/){
        $is_bcf = 1;
}


my $IN;
my $headers;
my $start_indiv_retrieval = 12;
my $chrom_retrieval = 2;
my $pos_retrieval = 3;
if ($is_vcf or $is_bcf){
	$start_indiv_retrieval = 9;
	$chrom_retrieval = 0;
	$pos_retrieval = 1;
	if ($is_vcf){
		$headers = `grep '#CHROM' $input`;
		open($IN,$input);
	}
	elsif ($is_bcf){
		my $cmd = "$VCFTOOLS_EXE --bcf $input --stdout --recode  | head -4000 | grep CHROM";
		$headers = `$cmd`;
		my $cmd2 = "$VCFTOOLS_EXE --bcf $input --stdout --recode";
		open $IN, '-|' , "$cmd2" or die "Can not run Vcftools";
	}
}
else{
	$headers= <$IN>;
	open($IN,$input);
}
$headers=~s/\n//g;
$headers=~s/\r//g;
my @ind_names = split(/\t/,$headers);
my @individual_names;
for (my $i = $start_indiv_retrieval; $i <= $#ind_names; $i++)
{
	push(@individual_names,$ind_names[$i]);
}
my %maximums;
while(<$IN>)
{
	my $line = $_;
	$line=~s/\n//g;
	$line=~s/\r//g;
	my @infos = split(/\t/,$line);
	if (scalar @infos > 8 && !/#CHROM/){
		my $chrom = $infos[$chrom_retrieval];
		my $position = $infos[$pos_retrieval];
		if ($position > $maximums{$chrom}){$maximums{$chrom}=$position;}
		my $classe_position = int($position/$step);
		$counts{$chrom}{$classe_position}++;
		
		my $ref_allele = $infos[11];
		if ($is_vcf or $is_bcf){
			$ref_allele = "0/0";
		}
		for (my $i = $start_indiv_retrieval; $i <= $#infos; $i++){
			if (!$counts_by_ind{$chrom}{$classe_position}{$i}){$counts_by_ind{$chrom}{$classe_position}{$i} = 0;}
			if ($infos[$i] ne $ref_allele){
				$counts_by_ind{$chrom}{$classe_position}{$i}++;
			}
		}
	}
}
close($IN);

#######################################################
# global
#######################################################
open(my $OUT,">$out");
print $OUT "Chromosome	Position	SNPs\n";
my $chr_num = 0;
foreach my $chrom(sort keys(%counts))
{
	$chr_num++;
	my $ref_counts = $counts{$chrom};
	my %final_counts = %$ref_counts;
	my $x = 0;
	#foreach my $classe_position(sort {$a<=>$b} keys(%final_counts))
	for (my $classe_position = 0; $classe_position <= $maximums{$chrom}/$step;$classe_position++)
	{
		my $nb = 0;
		if ($counts{$chrom}{$classe_position})
		{
			$nb = $counts{$chrom}{$classe_position};
		}
		$x += $step;
		print $OUT "$chrom	$x	$nb\n";
	}
	if ($chr_num >= $max_chr_num){last;}
}
close($OUT);

#######################################################
# For each individual
#######################################################
open(my $OUT2,">$out.by_sample");
$chr_num = 0;
print $OUT2 "Chromosome	".join("\t",@individual_names) . "\n";
foreach my $chrom(sort keys(%counts_by_ind))
{
	$chr_num++;
        my $ref_counts = $counts_by_ind{$chrom};
        my %final_counts = %$ref_counts;
	for (my $classe_position = 0; $classe_position <= $maximums{$chrom}/$step;$classe_position++)
        {
			print $OUT2 "$chrom";
			my $num_ind = $start_indiv_retrieval;
			foreach my $indiv(@individual_names)
			{
				my $val = 0;
				
				if ($counts_by_ind{$chrom}{$classe_position}{$num_ind})
				{
					$val = $counts_by_ind{$chrom}{$classe_position}{$num_ind};
				}
				print $OUT2 "	$val";
				$num_ind++;
			}
			print $OUT2 "\n";
        }
	if ($chr_num >= $max_chr_num){last;}
}
close($OUT2);
