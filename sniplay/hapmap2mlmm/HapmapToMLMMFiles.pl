#!/usr/bin/perl

use strict;
use Switch;
use Getopt::Long;

my $usage = qq~Usage:$0 <args> [<opts>]
where <args> are:
    -h, --hapmap       <Hapmap input file>
    -m, --map          <Map output file>
    -g, --geno         <Genotype output file>
    -p, --path         <Path for transpose executable>
~;
$usage .= "\n";

my ($hapmap,$map,$geno,$path);


GetOptions(
        "geno=s"       => \$geno,
        "map=s"        => \$map,
        "hapmap=s"     => \$hapmap,
        "path=s"       => \$path,
);


die $usage
  if ( !$geno || !$map || !$hapmap || !$path);
  
my $TRANSPOSE_EXE = "$path/transpose.awk";

my @snps;
my %chrom_pos;
my $num_line = 0;
open(my $O,">geno_transposed");
open(my $H,$hapmap);
while(<$H>)
{
	$num_line++;
	my $line = $_;
	chomp($line);
	$line =~s/\r//g;
	$line =~s/\n//g;
	my @infos = split(/\t/,$line);
	if ($num_line == 1)
	{
		print $O "Ind_id";
		for (my $i = 11; $i <= $#infos; $i++)
		{
			my $individual = $infos[$i];
			print $O "	" . $individual;
		}
		print $O "\n";
	}
	elsif ($num_line > 1)
	{
		my $snp = $infos[0];
		my $variation = $infos[1];
		my %scores;
		if ($variation =~/(\w)\/(\w)/)
		{
			my $allele1 = $1;
			my $allele2 = $2;
			$scores{$allele1} = 0;
			$scores{$allele2} = 1;
		}
		my $chrom = $infos[2];
		my $pos = $infos[3];
		$chrom_pos{$snp}{"chrom"} = $chrom;
		$chrom_pos{$snp}{"pos"} = $pos;
		push(@snps,$snp);
		print $O "$snp";
		for (my $i = 11; $i <= $#infos; $i++)
		{
			my $genotype = $infos[$i];
			my @alleles = split("",$genotype);
			if ($genotype ne "NN")
			{
				my $score = $scores{$alleles[0]} + $scores{$alleles[1]};
				print $O "	$score";
			}
			else
			{
				print $O "	NA";
			}
		}
		print $O "\n";
	}
}
close($H);
close($O);

open(my $M,">$map");
print $M "SNP	Chr	Pos\n";
foreach my $snp(@snps)
{
	print $M "$snp	" . $chrom_pos{$snp}{"chrom"} . "	". $chrom_pos{$snp}{"pos"} . "\n";
}
close($M);

system("$TRANSPOSE_EXE geno_transposed >geno_transposed2");

open(my $F,">$geno");
open(my $G,"geno_transposed2");
while(<$G>)
{
	my $line = $_;
	$line =~s/ /\t/g;
	print $F $line;
}
close($G);
close($F);

unlink("geno_transposed");
unlink("geno_transposed2");


