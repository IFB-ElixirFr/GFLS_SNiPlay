
#!/usr/bin/perl

use strict;
use Getopt::Long;
use Bio::SeqIO;

my $usage = qq~Usage:$0 <args> [<opts>]

where <args> are:

    -i, --input          <VCF input>
    -o, --out            <output basename>
~;
$usage .= "\n";

my ($input,$out);

GetOptions(
	"input=s"        => \$input,
	"out=s"          => \$out
);


die $usage
  if ( !$input);



my $nb_gene = `grep -c mRNA $input`;
$nb_gene =~s/\n//g;
my $nb_intergenic = `grep -c INTERGENIC $input`;
$nb_intergenic =~s/\n//g;

my $nb_intron = `grep -c INTRON $input`;
$nb_intron =~s/\n//g;
my $nb_UTR = `grep -c UTR $input`;
$nb_UTR =~s/\n//g;
my $nb_exon = $nb_gene - $nb_intron - $nb_UTR;

my $nb_ns = `grep -c NON_SYNONYMOUS_CODING $input`;
$nb_ns =~s/\n//g;
my $nb_s = $nb_exon - $nb_ns;



	
#system("$VCFTOOLS_EXE --vcf $input --remove-filtered-all --out $out --hardy >>vcftools.log 2>&1");
system("vcftools --vcf $input --remove-filtered-all --out $out --het >>vcftools.log 2>&1");
system("vcftools --vcf $input --remove-filtered-all --out $out --TsTv-summary >>vcftools.log 2>&1");
system("vcftools --vcf $input --remove-filtered-all --out $out --missing-indv >>vcftools.log 2>&1");

open(my $G,">$out.annotation");
print $G "Genic	$nb_gene\n";
print $G "Intergenic	$nb_intergenic\n";
print $G "========\n";
print $G "Intron	$nb_intron\n";
print $G "Exon	$nb_exon\n";
print $G "UTR	$nb_UTR\n";
print $G "========\n";
print $G "Non-syn	$nb_ns\n";
print $G "Synonym	$nb_s\n";
close($G);





	
	
