#!/usr/bin/perl

use strict;

my $vcf = $ARGV[0];

my %chrs;
my $ok = 0;
open(my $V,$vcf);
while(<$V>)
{
	if ($ok)
	{
		my ($chr,$pos) = split(/\t/,$_);
		$chrs{$chr}++;
	}
	if (/#CHROM/){$ok = 1;}
}
close($V);

foreach my $chr(sort keys(%chrs))
{
	my $nb = $chrs{$chr};
	print "$chr	$nb\n";
}
