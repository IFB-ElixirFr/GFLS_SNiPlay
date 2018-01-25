
#!/usr/bin/perl

use strict;
use Getopt::Long;

my $usage = qq~Usage:$0 <args> [<opts>]

where <args> are:

    -i, --in           <PED input>
    -o, --out          <Fasta output>
~;
$usage .= "\n";

my ($input,$out);



GetOptions(
	"in=s"          => \$input,
	"out=s"          => \$out,
);


die $usage
  if ( !$input || !$out);
  

my %IUPAC =
(
	'00'=> "?",
	'AA'=> "A",
	'CC'=> "C",
	'GG'=> "G",
	'TT'=> "T",
        'AG'=> "R",
        'GA'=> "R",
        'CT'=> "Y",
        'TC'=> "Y",
        'TG'=> "K",
        'GT'=> "K",
        'CG'=> "S",
        'GC'=> "S",
        'AT'=> "W",
        'TA'=> "W",
        'AC'=> "M",
        'CA'=> "M",
);

open(my $O,">$out");
open(my $P,$input) or die "File does not exist";
while(<$P>)
{
	my $line = $_;
	$line =~s/\r//g;
	$line =~s/\n//g;
	my @infos = split("\t",$_);
	my $ind = $infos[0];
	print $O ">$ind\n";
	for (my $i = 6; $i <= $#infos; $i= $i+2)
	{
		my $code = $infos[$i].$infos[$i+1];
		my $letter = $IUPAC{$code};
		print $O $letter;
	}
	print $O "\n";
}
close($P);
close($O);
