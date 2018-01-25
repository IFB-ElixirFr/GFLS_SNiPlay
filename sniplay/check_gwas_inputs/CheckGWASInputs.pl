#!/usr/bin/perl

use strict;
use Switch;
use Getopt::Long;

my $usage = qq~Usage:$0 <args> [<opts>]
where <args> are:
    -h, --hapmap       <Hapmap input file>
    -t, --trait        <Trait input file>
    -o, --out          <Output base name>
~;
$usage .= "\n";

my ($hapmap,$trait,$out);


GetOptions(
        "trait=s"      => \$trait,
        "out=s"        => \$out,
        "hapmap=s"     => \$hapmap
);


die $usage
  if ( !$trait || !$out || !$hapmap);
  
my %inds;

#######################################
# get individuals in trait file
#######################################
my %traits;
my $head_trait = `head -1 $trait`;
open(my $T,$trait);
<$T>;
while(<$T>)
{
	my @infos = split(/\t/,$_);
	my $ind = $infos[0];
	$inds{$ind}++;
	$traits{$ind} = $_;
}
close($T);
my $nb_ind_trait = scalar keys(%traits);

#######################################
# get individuals in hapmap file
#######################################
my $line_ind = `head -1 $hapmap`;
chomp($line_ind);
my @infos = split(/\t/,$line_ind);
for (my $i = 11; $i <= $#infos; $i++)
{
	my $ind = $infos[$i];
	$inds{$ind}++;
}
my $nb_ind_hapmap = scalar @infos - 11;

#################################################################
# create trait output by keeping individuals found in both files
#################################################################
open(my $O,">$out.trait");
print $O $head_trait;
my $nb_common = 0;
foreach my $ind(keys(%inds))
{
	my $nb_found = $inds{$ind};
	if ($nb_found == 2)
	{
		$nb_common++;
		print $O $traits{$ind};
	}
}
close($O);


#####################################################################
# create hapmap output after keeping individuals found in both files
# and removing monomorphic positions
#####################################################################
open(my $O2,">$out.hapmap");
my $numline = 0;
my %genotypes;
my %columns_to_keep;
my $nb_monomorphic = 0;
my $not_biallelic = 0;
my $diff_variation = 0;
open(my $H,$hapmap);
while(<$H>)
{
	$numline++;
	my $line = $_;
	$line =~s/\n//g;
	$line =~s/\r//g;
	my @infos = split(/\t/,$line);
	if ($numline == 1)
	{
		my @titles;
		for (my $i = 0; $i <= 10; $i++)
		{
			my $title = $infos[$i];
			push(@titles,$title);
		}
		print $O2 join("\t",@titles);
		for (my $i = 11; $i <= $#infos; $i++)
		{
			my $ind = $infos[$i];
			my $nb_found = $inds{$ind};
			if ($nb_found == 2)
			{
				print $O2 "	$ind";
				$columns_to_keep{$i} = 1;
			}
		}
		print $O2 "\n";
	}
	else
	{
		my $to_be_printed = "";
		my $variation = $infos[1];
        for (my $i = 0; $i <= 10; $i++)
		{
			my $title = $infos[$i];
			$to_be_printed .= "$title	";
		}
		my %letters;
		for (my $i = 11; $i <= $#infos; $i++)
		{
			if ($columns_to_keep{$i})
			{
				my $genotype = $infos[$i];
				if ($genotype ne 'NN')
				{
					my ($allele1,$allele2) = split(//,$genotype);
					$letters{$allele1}=1;
					$letters{$allele2}=1;
				}
				$to_be_printed .= "$genotype	";
			}
		}
		chop($to_be_printed);
		
		my $variation_obs = join("/",sort keys(%letters));
		
		# print only if polymorphic
		if (scalar keys(%letters) < 2)
		{
			$nb_monomorphic++;
		}
		elsif (scalar keys(%letters) > 2)
		{
			$not_biallelic++;
		}
		else
		{
			if ($variation ne $variation_obs)
			{
				$to_be_printed =~s/$variation/$variation_obs/;
				$diff_variation++;
			}
			
			print $O2 $to_be_printed . "\n";
		}
	}
}
close($H);
close($O2);

print "==============================================\n";
print "Individuals\n";
print "==============================================\n";
print "Individuals in hapmap file: $nb_ind_hapmap\n";
print "Individuals in trait file: $nb_ind_trait\n";
print "Individuals found in both files: $nb_common\n";
print "==============================================\n";
print "Markers\n";
print "==============================================\n";
print "Discarded markers:\n";
print "Monomorphic: $nb_monomorphic\n";
print "Not biallelic: $not_biallelic\n";
print "Modified markers:\n";
print "Difference in variation: $diff_variation\n";

