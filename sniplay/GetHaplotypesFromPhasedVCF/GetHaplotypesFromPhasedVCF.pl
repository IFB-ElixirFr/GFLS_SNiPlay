#!/usr/bin/perl

use strict;

my $vcf = $ARGV[0];
my $out = $ARGV[1];

open(O1,">$out.haplo.fas");
open(O2,">$out.distinct_haplotypes.txt");
open(O3,">$out.distinct_haplotypes.fa");

my %indiv;
my %genes;
my $nb_cols = 0;
my %phasing;
open(my $V,$vcf);
while(<$V>){
	my $line = $_;
	$line =~s/\n//g;
	$line =~s/\r//g;
	my @infos = split(/\t/,$line);
	if (/#CHROM/){
		for (my $i = 9; $i <= $#infos; $i++){
			$indiv{$i} = $infos[$i];		
		} 
	}
	else{
		my $gene = $infos[0];
		my $ref = $infos[3];
		my $alt = $infos[4];
		$nb_cols = $#infos;
		for (my $i = 9; $i <= $#infos; $i++){
			my ($alleles,$depth) = split(":",$infos[$i]);
			$alleles =~s/0/$ref/g;
			$alleles =~s/1/$alt/g;
			my $ind = $indiv{$i};
			#if ($alleles =~/\// && $genes{$gene}){next;}
			$phasing{$gene}{$i}.= $alleles ." ";
		}
		$genes{$gene} = 1;
	}
}
close($V);

my %haplos;
my $gene_ok = 0;
my $gene_shared = 0;
my $gene_all_shared = 0;
my $nb_gene = 0;
my %haplotypes2;
foreach my $gene(keys(%phasing))
{
	my $list = "";
	my $ok = 0;
	my %haplotypes;
	for (my $i=9;$i <= $nb_cols;$i++){
		
		my $alleles = $phasing{$gene}{$i};
		if ($alleles eq "" or $alleles =~ /\./ or $alleles =~/ \w+\//){next;}
		my @snps = split(" ",$alleles);
		if (scalar @snps < 2){next;}
		$alleles =~s/\//\|/g;	
		my @al = split(" ",$alleles);
		my $haplo1 = "";
		my $haplo2 = "";
		foreach my $a(@al){
			my ($a1,$a2) = split(/\|/,$a);
			$haplo1.= $a1;
			$haplo2.= $a2;
		}
		$haplotypes{$haplo1}++;
		$haplotypes{$haplo2}++;
		my $ind = $indiv{$i};
		$haplotypes2{$gene}{$haplo1}.=$ind."_1,";
		$haplotypes2{$gene}{$haplo2}.=$ind."_2,";
		$haplos{$gene}{$haplo1}++;
		$haplos{$gene}{$haplo2}++;
		$list .= ">".$gene."_".$ind."_1\n$haplo1\n";
		$list .= ">".$gene."_".$ind."_2\n$haplo2\n";
		$ok++;
	}
	#print "$gene $nb_cols $ok\n";
	#if ($ok > 13 && $found_M1_M2 == 2){
	if ($ok == ($nb_cols - 8)){
		$nb_gene++;
		print O1 $list;
	}
}
foreach my $gene(sort {$a<=>$b} keys(%haplos)){
	print O2 "===$gene===\n";
	my $ref_hash = $haplos{$gene};
	my %hash = %$ref_hash;
	my $num_haplo = 0;
	foreach my $haplo(keys(%hash)){
		$num_haplo++;
		my $haplo_name = "haplo".$num_haplo;
		my $nb = $haplos{$gene}{$haplo};
		my $ind = $haplotypes2{$gene}{$haplo};
		print O2 $haplo_name.":$nb:".$haplotypes2{$gene}{$haplo}."\n".$haplo."\n"; 
		if ($nb >= 1){
			print O3 ">".$haplo_name."|$nb\n";
                        print O3 $haplo."\n";
		}
	}
}

close(O1);
close(O2);
close(O3);
#print scalar keys(%haplos);

