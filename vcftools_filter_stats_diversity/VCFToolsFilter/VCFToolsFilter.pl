
#!/usr/bin/perl

use strict;
use Getopt::Long;
use Bio::SeqIO;

my $usage = qq~Usage:$0 <args> [<opts>]

where <args> are:

    -i, --input          <VCF input>
    -o, --out            <Output basename>
      
      <opts> are:

    -s, --samples        <Samples to be analyzed. Comma separated list>
    -c, --chromosomes    <Chromosomes to be analyzed. Comma separated list>
    -e, --export         <Output format (VCF/freq/plink. Default: VCF>
    -f, --frequency      <Minimum MAF. Default: 0.001>
    -m, --max_freq       <Maximum MAF. Default: 0.5>
    -a, --allow_missing  <Allowed missing data proportion per site. Must be comprised between 0 and 1. Default: 1>
    -n, --nb_alleles     <Accepted number of alleles (min,max). Default: 2,4>
    -t, --type           <Type of polymorphisms to keep (ALL/SNP/INDEL). Default: ALL>
    -b, --bounds         <Lower bound and upper bound for a range of sites to be processed (start,end). Default: 1, 100000000>
~;
$usage .= "\n";

my ($input,$out);


#my $indel_size_max = 500;
#my $indel_size_min = 1;
my $frequency_max = 0.5;
my $frequency_min = 0.001;
my $pos_max = 100000000000;
my $pos_min = 0;
my $filter_snp_type = "all";

my $missing_data = 1;
my $export = "VCF";
my $type = "ALL";
my $nb_alleles;
my $bounds;
my $samples;
my $chromosomes;

GetOptions(
	"input=s"        => \$input,
	"out=s"          => \$out,
	"samples=s"      => \$samples,
	"chromosomes=s"  => \$chromosomes,
	"frequency=s"    => \$frequency_min,
	"max_freq=s"     => \$frequency_max,
	"allow_missing=s"=> \$missing_data,
	"export=s"       => \$export,
	"type=s"         => \$type,
	"nb_alleles=s"   => \$nb_alleles,
	"bounds=s"       => \$bounds,
);


die $usage
  if ( !$input || !$out);


my @dnasamples;
if ($samples)
{
	@dnasamples = split(",",$samples);
}
my @nalleles;
if ($nb_alleles)
{
	@nalleles = split(",",$nb_alleles);
}
my @boundaries;
if ($bounds)
{
	@boundaries = split(",",$bounds);
}
my @chromosomes_list;
if ($chromosomes)
{
	@chromosomes_list = split(",",$chromosomes);
}


my $experiment = "chromosomes";
my $table = "";
my %genes;
my @snp_ids;
my @snp_ids_and_positions;
my @snp_ids_and_positions_all;
my $gene;
my $snp_num = 0;
my %ref_sequences;
my %snps_of_gene;

	


my $indiv_cmd = "";
if (@dnasamples)
{
	$indiv_cmd = "--indv " . join(" --indv ",@dnasamples);
}

my $chrom_cmd = "";
if (@chromosomes_list)
{
	$chrom_cmd = "--chr " . join(" --chr ",@chromosomes_list);
}

my $export_cmd = "--recode";
if ($export eq "freq")
{
	$export_cmd = "--freq";
}
if ($export eq "plink")
{
	$export_cmd = "--plink";
}
 


my $nb_alleles_cmd = "--min-alleles 1 --max-alleles 4";
if (@nalleles)
{
	$nb_alleles_cmd = "--min-alleles $nalleles[0] --max-alleles $nalleles[1]";
}
my $bounds_cmd = "--from-bp 1 --to-bp 100000000";
if (@boundaries)
{
        $bounds_cmd = "--from-bp $boundaries[0] --to-bp $boundaries[1]";
}

 
my $type_cmd = "";
if ($type eq "INDEL")
{
	$type_cmd = "--keep-only-indels";
}
if ($type eq "SNP")
{
	$type_cmd = "--remove-indels";
}

	
system("vcftools --vcf $input --out $out --keep-INFO-all --remove-filtered-all $type_cmd $export_cmd $chrom_cmd $indiv_cmd $nb_alleles_cmd --maf $frequency_min --max-maf $frequency_max --max-missing $missing_data >>vcftools.log 2>&1");





	
	
