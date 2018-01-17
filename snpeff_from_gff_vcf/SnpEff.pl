#!/usr/bin/perl

use strict;
use Getopt::Long;

my $usage = qq~Usage:$0 <args> [<opts>]
where <args> are:
    -i, --input         <input VCF>
    -o, --output        <output>
    -g, --gff           <GFF annotation>
    -f, --fasta         <Fasta of chromosomes>
    -h, --html          <HTML output>
~;
$usage .= "\n";

my ($input,$output,$gff,$fasta,$html);


GetOptions(
	"input=s"      => \$input,
	"output=s"     => \$output,
	"gff=s"        => \$gff,
	"fasta=s"      => \$fasta,
	"html=s"       => \$html
);


die $usage
  if ( !$input || !$output || !$fasta || !$gff || !$html);


if (!-e $gff){
        die "Error: GFF input does not exist\n"
}
if (!-e $fasta){
        die "Error: Fasta input does not exist\n"
}

#my $SNPEFF_PATH = "/usr/local/bioinfo/galaxy/galaxy_dist/tools/SNiPlay/SnpEff/snpEff";
my $SNPEFF_PATH = $ENV{SNPEFF_JAR_PATH};


my $session = $$;
mkdir($session);
mkdir("$session/data");
mkdir("$session/data/genomes");
mkdir("$session/data/myspecies");

system("cp -rf $fasta $session/data/genomes/myspecies.fa");
system("cp -rf $gff $session/data/myspecies/genes.gff");

open(my $C,"$SNPEFF_PATH/snpEff.config");
open(my $C2,">$session/snpEff.config");
while(<$C>)
{
	if (/data_dir/)
	{
		print $C2 "data_dir = ./data\n";
	}
	elsif (/^genomes/)
	{
		print $C2 "genomes : \\n";
        	print $C2 "myspecies, myspecies \\n";
	}
	else
	{
		print $C2 $_;
	}
}
print $C2 "myspecies.genome : myspecies\n";
close($C);
close($C2);


my $build_cmd = "snpEff build -c $session/snpEff.config -gff3 myspecies";
#my $build_cmd = "java -jar $SNPEFF_PATH/snpEff.jar build -c $session/snpEff.config -gff3 myspecies";
system($build_cmd);

my $eff_cmd = "snpEff eff -c $session/snpEff.config -o vcf -no-downstream -no-upstream myspecies -s $html $input >$output";
#my $eff_cmd = "java -jar $SNPEFF_PATH/snpEff.jar eff -c $session/snpEff.config -o vcf -no-downstream -no-upstream myspecies -s $html $input >$output";
system($eff_cmd);


system("rm -rf $session");
