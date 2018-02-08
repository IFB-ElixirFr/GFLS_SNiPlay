#!/usr/bin/perl

use strict;
use Getopt::Long;
use Bio::SeqIO;
use Cwd ;
use FindBin qw ( $Bin $Script );

my $CURRENT_DIR = $Bin;

my $ROOTING_EXE = "java -jar ". $CURRENT_DIR . "/Rootings_54.jar";

my $usage = qq~Usage:$0 <args> [<opts>]
where <args> are:
    -i, --input         <newick input>
    -o, --output        <newick output>
~;
$usage .= "\n";

my ($input,$outfile);


GetOptions(
	"input=s"    => \$input,
	"output=s"   => \$outfile
);


die $usage
  if ( !$input || !$outfile);
  
my $treefile = $input;
	

# replace negative values by 0
open(T,$treefile);
open(T2,">$treefile.2");
while(<T>)
{
	my $line = $_;
	$line =~s/\-\d+\.*\d*\,/0,/g;
	$line =~s/\-\d+\.*\d*\)/0\)/g;
	print T2 $line;
}
close(T);
close(T2);
	
my $rooting_command = $ROOTING_EXE . " -input $treefile.2 -output $treefile.all -midpoint $treefile.midpoint >>$treefile.rooting.log 2>&1";
system($rooting_command);

unlink("$treefile.all");
unlink("$treefile.2");
rename("$treefile.midpoint",$outfile);

	
	


