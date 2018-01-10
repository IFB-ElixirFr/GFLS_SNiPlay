
#!/usr/bin/perl

use strict;
use Getopt::Long;
use Bio::SeqIO;

my $usage = qq~Usage:$0 <args> [<opts>]

where <args> are:

    -i, --input          <VCF input>
    -o, --out            <output basename>
    -w, --window         <window size>

<opts> are:
    -g, --group          <group file>
~;
$usage .= "\n";

my ($input,$out);
my $window = 200000;
my $groupfile;
GetOptions(
	"input=s"        => \$input,
	"out=s"          => \$out,
	"window=s"       => \$window,
	"group=s"        => \$groupfile
);



die $usage
  if ( !$input);

my %hash;
my $cmd_part = "";
if ($groupfile && -e $groupfile)
{
	open(my $G,$groupfile) or die "Cannot open $groupfile: $!";
	while(<$G>)
	{
		my $line = $_;
		chomp($line);
		$line=~s/\r//g;
		$line=~s/\n//g;
		my @infos = split(/;/,$line);
		if ($infos[0] && $infos[1])
		{
			$hash{$infos[1]} .= " --indv " . $infos[0];
			$cmd_part .= " --indv " . $infos[0];
		}
	}
	close($G);
}


if ($window =~/^(\d+)\s*$/){
	$window = $1;
}
else{
	die "Error: window size must be an integer\n";
}


my $input_format = "--vcf $input.sorted";
if ($input =~/\.vcf/){
	system("vcf-sort $input >$input.sorted");
	$input_format = "--vcf $input.sorted";
}
elsif ($input =~/\.bcf/){
	$input_format = "--bcf $input";
}
system("vcftools $input_format --out $out --window-pi $window --window-pi-step $window >>$out.vcftools.log 2>&1");
system("vcftools $input_format --out $out --TajimaD $window >>$out.vcftools.log 2>&1");
system("vcftools $input_format --out $out --TsTv $window >>$out.vcftools.log 2>&1");
system("vcftools $input_format --out $out --SNPdensity $window >>$out.vcftools.log 2>&1");



##################################################
# Tajima distribution
##################################################
my %counts3;
my @values;
open(my $TAJIMA,"$out.Tajima.D");
<$TAJIMA>;
my $max = 0;
while(<$TAJIMA>){
        my $line = $_;
        $line =~s/\n//g;
        $line =~s/\r//g;
        my @inf = split(/\t/,$line);
        my $val = $inf[3];
	my $arrondi = sprintf("%.1f",$val);
        $counts3{$arrondi}++;
}
close($TAJIMA);


my %counts3_pop;

if (keys(%hash) > 0)
{
	my $files_pi = "";
	my $files_dtajima = "";
	my $cmd_fst = "--fst-window-size $window";
	my $cmd_fst_by_marker = "";
	foreach my $pop(sort(keys(%hash)))
	{
		my $cmd_part = $hash{$pop};

		
		open(my $POP,">$input.$pop.txt");
		my @ind_of_pop = split(" --indv ",$cmd_part);
		print $POP join("\n",@ind_of_pop);
		close($POP);
		$cmd_fst .= " --weir-fst-pop $input.$pop.txt";
		$cmd_fst_by_marker .= " --weir-fst-pop $input.$pop.txt";

        	my $cmd_part = $hash{$pop};
		system("vcftools $input_format --remove-filtered-all --out $out.$pop --window-pi $window --window-pi-step $window $cmd_part --maf 0.001 >>$out.vcftools.log 2>&1");
		my $sed_cmd = "sed -i \"s\/PI\/$pop\/g\" $out.$pop.windowed.pi";
		system($sed_cmd);
		$files_pi .= "$out.$pop.windowed.pi ";

		system("vcftools $input_format --remove-filtered-all --out $out.$pop --SNPdensity $window $cmd_part --maf 0.001 >>$out.vcftools.log 2>&1");

		system("vcftools $input_format --remove-filtered-all --out $out.$pop --TajimaD $window $cmd_part --maf 0.001 >>$out.vcftools.log 2>&1");
		my $sed_cmd = "sed -i \"s\/TajimaD\/$pop\/g\" $out.$pop.Tajima.D";
                system($sed_cmd);
		$sed_cmd = "sed -i \"s/nan/0/g\" $out.Tajima.D";
                system($sed_cmd);
                $files_dtajima .= "$out.$pop.Tajima.D ";
		

		open(my $TAJIMAPOP,"$out.$pop.Tajima.D");
		<$TAJIMAPOP>;
		my $max = 0;
		while(<$TAJIMAPOP>){
			my $line = $_;
			$line =~s/\n//g;
			$line =~s/\r//g;
			my @inf = split(/\t/,$line);
			my $val = $inf[3];
			my $arrondi = sprintf("%.1f",$val);
			$counts3_pop{$pop}{$arrondi}++;
		}
		close($TAJIMAPOP);
	
		system("vcftools $input_format --remove-filtered-all --out $out.$pop --TsTv $window $cmd_part --maf 0.001 >>$out.vcftools.log 2>&1");
	}
	system("paste $files_pi >$out.combined.pi");
	my $awk_cmd = "awk {'print \$1\"\t\"\$2\"\t\"\$5\"\t\"\$10'} $out.combined.pi >$out.combined.pi.txt";
	system($awk_cmd);

	system("paste $files_dtajima >>$out.combined.dtajima");
        $awk_cmd = "awk {'print \$1\"\t\"\$2\"\t\"\$4\"\t\"\$8'} $out.combined.dtajima >$out.combined.dtajima.txt";
        system($awk_cmd);

	system("vcftools $input_format --remove-filtered-all --out $out.fst $cmd_fst  --maf 0.001 >>$out.vcftools.log 2>&1");
	$awk_cmd = "awk {'print \$1\"\t\"\$2\"\t\"\$5'} $out.fst.windowed.weir.fst >$out.fst.txt";
        system($awk_cmd);

	system("vcftools $input_format --remove-filtered-all --out $out.fst $cmd_fst_by_marker --maf 0.001 >>$out.vcftools.log 2>&1");

	my %genes;
	open(my $V,$input);
	my $go = 0;
	while(<$V>){
		if ($go){
			my @infos = split(/\t/,$_);
			my $chr = $infos[0];
			my $pos = $infos[1];
			my $annot = $infos[7];
			my $gene = "#";
			if ($annot =~/EFF=\w+\((.*)\)/){
				my @parts = split(/\|/,$1);
				$gene = $parts[4];
			}
			$genes{"$chr:$pos"} = $gene;
		}
		if (/^#CHROM/){$go = 1;}
	}

	my $cumul = 0;
	my $previous_pos;
	my %chrs;
	
	# get only the first 2000 points
	`sort -k 3 $out.fst.weir.fst | grep -v 'nan' | grep -v 'e-' | tail -10001 | sort -n -k 2 >$out.fst.weir.highest.fst`;
	
	my %hash;
	open(my $F,"$out.fst.weir.highest.fst");
	<$F>;
	open(my $O,">$out.fst.by_marker.txt");
	open(my $O2,">$out.fst.by_marker.genes.txt");
	print $O "Chromosome	Marker	Position	Fst\n";
	print $O2 "Chromosome	Marker	Position	Fst\n";
	while(<$F>){
		my $line = $_;
		chomp($line);
		my ($chr,$pos,$val) = split(/\t/,$line);
		$hash{$chr}{$pos} = $val;
	}
	close($F);
	foreach my $chr(sort keys(%hash)){
		my $ref_hash = $hash{$chr};
		my %hash2 = %$ref_hash;
		foreach my $pos(sort {$a<=>$b} keys(%hash2)){
			my $val = $hash{$chr}{$pos};
			if (!$chrs{$chr}){
				$cumul = $previous_pos;
				$chrs{$chr} = 1;
			}
			my $x = $pos+$cumul;
			#print $O2 "$chr	$chr:$pos	$pos	$val	" . $genes{"$chr:$pos"} . "\n";
			print $O2 "$chr	$chr:$pos	$pos	$val\n";
			print $O "$chr	$chr:$pos	$x	$val\n";
			$previous_pos = $x;
		}
		
	}
	
	close($O);
	close($O2);

	#$awk_cmd = "awk {'if (\$3>0.6 && \$3<=1 && \$3!~\"nan\"){print \$1\"\t\"\$2\"\t\"\$2\"\t\"\$3}'} $out.fst.weir.fst >$out.fst.by_marker.txt";
        #system($awk_cmd);
}
open(my $M3,">$out.tajima_distrib");
print $M3 "Values	nb of values";
if (keys(%hash) > 0){
	print $M3 "	".join("\t",keys(%hash));
}
print $M3 "\n";
for (my $i = -2.5; $i <= 2.5; $i=$i+0.1)
{
        my $arrondi = sprintf("%.1f",$i);
        my $nb = 0;
        if ($counts3{$arrondi})
        {
                $nb = $counts3{$arrondi};
        }
	print $M3 "Tajima	$arrondi	$nb";
	if (keys(%hash) > 0){
		foreach my $pop(sort(keys(%hash))){
			my $nb_pop = 0;
			if ($counts3_pop{$pop}{$arrondi}){
				$nb_pop = $counts3_pop{$pop}{$arrondi};
			}
			print $M3 "	$nb_pop";
		}
	}
        print $M3 "\n";
}
close($M3);



	
	
