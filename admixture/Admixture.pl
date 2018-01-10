#!/usr/bin/perl

use strict;
use Getopt::Long;
use File::Basename;

my $usage = qq~Usage:$0 <args> [<opts>]
where <args> are:
    -i, --input         <input basename (for PED, FAM, BIM)>
    -o, --output        <output>
    -k, --kmin          <K min. int>
    -m, --maxK          <K max. int>
    -d, --directory     <temporary directory>
    -t, --threshold     <threshold admixture proportion for group assignation>
~;
$usage .= "\n";

my ($input,$output,$kmin,$kmax,$directory,$threshold);


GetOptions(
	"input=s"      => \$input,
	"output=s"     => \$output,
	"kmin=s"       => \$kmin,
	"maxK=s"       => \$kmax,
	"directory=s"  => \$directory,
	"threshold"    => \$threshold,
);


die $usage
  if ( !$input || !$output || !$kmin || !$kmax || !$directory );

if ($kmin =~/^(\d+)\s*$/){
        $kmin = $1;
}
else{
        die "Error: kmin must be an integer\n";
}
if ($kmax =~/^(\d+)\s*$/){
        $kmax = $1;
}
else{
        die "Error: kmax must be an integer\n";
}



open(my $P2,"$input.fam");
my $n = 0;
my $ind_num = 0;
my @individus;
while(<$P2>)
{
	my ($ind,$other) = split(/ /,$_);
	push(@individus,$ind);
}
close($P2);

my $basename = basename($input);

###################################
# launch admixture for different K
###################################
my %errors;
for (my $k = $kmin; $k <= $kmax; $k++)
{
	system("admixture --cv $input.bed $k >>$directory/log.$k 2>&1");
	my $cv_error_line = `grep -h CV $directory/log.$k`;
	if ($cv_error_line =~/: (\d+\.*\d*)$/)
	{
		$errors{$1} = $k;
	}
	system("cat $directory/log.$k >>$directory/logs");
	system("echo '\n\n====================================\n\n' >>$directory/logs");

	open(my $O2,">$basename.$k.final.Q");
	open(my $O3,">$directory/groups.$k");
	open(my $O,"$basename.$k.Q");
	my %hash_groupes;
	my %hash_indv;
	my %group_of_ind;
	my $i = 0;
	while (<$O>){
		$i++;
		my $line = $_;
		$line =~s/\n//g;
		$line =~s/\r//g;
		my @infos = split(/\s+/,$line);
		my $group = "admix";
		my $ind = $individus[$i];
		for (my $j = 0; $j <$k; $j++){
			my $val = $infos[$j];
			if ($val > ($threshold/100)){$group = "Q$j";}
		}
		if ($ind){      
			$hash_indv{$ind} = join("	",@infos);
			$hash_groupes{$group}{"ind"} .= ",".$ind;
			$group_of_ind{$ind} = $group;
		}
	}
	close($O);

	foreach my $group(sort keys(%hash_groupes)){
		my @inds = split(",",$hash_groupes{$group}{"ind"});
		foreach my $ind(@inds){
			if ($ind =~/\w+/){
				print $O3 "$ind;$group\n";
				print $O2 $ind."	".$hash_indv{$ind}. "\n";
			}
		}
	}
	close($O3);
	close($O2);

	system("cat $basename.$k.final.Q >>$directory/outputs.Q");
	system("echo '\n\n====================================\n\n' >>$directory/outputs.Q");
	system("cat $basename.$k.P >>$directory/outputs.P");
	system("echo '\n\n====================================\n\n' >>$directory/outputs.P");
}

my @sorted_errors = sort {$a<=>$b} keys(%errors);
my $best_K = $errors{@sorted_errors[0]};


open(BEST1,"$basename.$best_K.final.Q");
open(BEST2,">$directory/output");
print BEST2 "<Covariate>\n";
print BEST2 "<Trait>";
for (my $j=1;$j<=$best_K;$j++)
{
	print BEST2 "	Q" . $j;
}
print BEST2 "\n";
my $i = 0;
while(<BEST1>)
{
	my $line = $_;
	$line =~s/ /\t/g;
	print BEST2 $line;
	$i++;
}
close(BEST1);
close(BEST2);

system("cp -rf $directory/log.$best_K $directory/log");
system("cp -rf $directory/groups.$best_K $directory/groups");



