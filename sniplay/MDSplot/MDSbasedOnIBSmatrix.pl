#!/usr/bin/perl

use strict;
use Switch;
use Getopt::Long;
use Bio::SeqIO;

my $PLINK_EXE= "plink";

my $usage = qq~Usage:$0 <args> [<opts>]
where <args> are:
    -i, --in         <input>
    -o, --out        <output>
~;
$usage .= "\n";

my ($in,$out);


GetOptions(
	"in=s"        => \$in,
	"out=s"       => \$out
);

die $usage
  if ( !$in || !$out);
  
	
my $plink_command = $PLINK_EXE . " --file $in --noweb --cluster --matrix --mds-plot 2 --out $out >>$in.plink.log 2>&1";
system($plink_command);

my $awk_cmd = "awk \{\'print \$1\'\} $in.ped";
my $inds = `$awk_cmd`;
my @individuals = split("\n",$inds);

my %populations;
if (-e "$in.individual_info.txt")
{
	open(my $I,"$in.individual_info.txt");
	while(<$I>)
	{
		my $line = $_;
		$line =~s/\n//g;
		$line =~s/\r//g;
		my ($ind,$pop) = split(/;/,$line);
		$populations{$ind} = $pop;
	}
	close($I);
}

open(my $OUT,">$out.mds_plot.txt");
my $go = 0;
open(my $O,"$out.mds");
while(<$O>)
{
	if ($go)
	{
		my $line = $_;
		$line =~s/\n//g;
		$line =~s/\r//g;
		my @i = split(/\s+/,$line);
		if ($line =~/^ /)
		{
			my $ind = $i[1];
			my $pop = "Pop1";
			if ($populations{$ind})
			{
				$pop = $populations{$ind};
			}
			print $OUT "$pop	$ind	".$i[4]."	".$i[5]."\n";
		}
		if ($line =~/^\w/)
		{
			my $ind = $i[0];
			my $pop = "Pop1";
			if ($populations{$ind})
			{
				$pop = $populations{$ind};
			}
			print $OUT "$pop	$ind	".$i[3]."	".$i[4]."\n";
		}
		
	}
	if (/C1/){$go = 1;}
}
close($O);
close($OUT);


my $j = 0;
open(my $IBS,">$out.ibs_matrix.txt");
print $IBS "Individuals	" . join("\t",@individuals)."\n";
open(my $O2,"$out.mibs");
while(<$O2>)
{
	my $line = $_;
	$line =~s/\n//g;
	$line =~s/\r//g;
	my @i = split(/\s+/,$line);
	print $IBS $individuals[$j]. "	". join("\t",@i)."\n";
	$j++;
}
close($O2);
close($IBS);


	
	


