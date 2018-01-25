#!/usr/bin/perl

use strict;
use Switch;
use Getopt::Long;
use Bio::SeqIO;


my $usage = qq~Usage:$0 <args> [<opts>]
where <args> are:
    -i, --input         <FASTA input>
    -o, --output        <output filename>
    -d, --directory     <directory of egglib package>
~;
$usage .= "\n";

my ($infile,$outfile,$dir_exe);


GetOptions(
	"input=s"    => \$infile,
	"output=s"   => \$outfile,
	"directory=s"=> \$dir_exe
);


die $usage
  if ( !$infile || !$outfile || !$dir_exe);
  

my $EGGSTATS_EXE = "$dir_exe/egglib-2.1.5/bin/eggstats";

my %gene_alignments;
my $in  = Bio::SeqIO->new(-file => $infile , '-format' => 'Fasta');
while ( my $seq = $in->next_seq() ) 
{
	my $id = $seq -> id();
	my $sequence = $seq -> seq();
	my ($gene,$ind,$num_allele) = split("_",$id);
	$gene_alignments{$gene}.= ">$id\n$sequence\n";
}

open(OUT,">$outfile");
foreach my $gene(keys(%gene_alignments))
{
	open(F,">$gene.egglib_input.fa");
	print F $gene_alignments{$gene};
	close(F);
	
	my $results_egglib = `$EGGSTATS_EXE $gene.egglib_input.fa`;

	# parse Seqlib output
	if ($results_egglib)
	{
		my %egglig_stats;
		my @eggstats = split(/^/,$results_egglib);
		foreach my $eggstat(@eggstats)
		{
			my ($desc,$value) = split(/: /,$eggstat);
			chomp($value);
			$egglig_stats{$desc} = $value;
		}
		print OUT "$gene;";
		print OUT $egglig_stats{"Total number of sequences"} . ";";
		print OUT $egglig_stats{"Total number of sites"} . ";";
		print OUT $egglig_stats{"Number of analyzed sites"} . ";";
		print OUT $egglig_stats{"S"} . ";";
		print OUT $egglig_stats{"thetaW"} . ";";
		print OUT $egglig_stats{"Pi"} . ";";
		print OUT $egglig_stats{"D"} . ";";
		print OUT $egglig_stats{"number of haplotypes"} . ";";
		print OUT $egglig_stats{"haplotypes diversity"} . ";";
		print OUT $egglig_stats{"Fay and Wu H"} . ";";
		print OUT $egglig_stats{"Fst"} . ";";
		print OUT $egglig_stats{"Snn"} . ";";
		print OUT "\n";
		unlink("$gene.egglib_input.fa");
	}
}
close(OUT);

