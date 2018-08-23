#!/usr/bin/perl

use strict;
use Getopt::Long;
use Bio::SeqIO;

use Cwd;
my $dir = getcwd;


my $usage = qq~Usage:$0 <args> [<opts>]
where <args> are:
    -i, --input         <input>
    -o, --output        <output>
    -d, --dotfile       <dotfile>
<opts> are:
    -g, --groups        <groupfile>
    -s, --stats         <statfile>
    -t, --tool_path     <tool_path>
~;
$usage .= "\n";

my ($infile,$output,$outfile,$groupfile,$statfile,$tool_path);


GetOptions(
	"input=s"    => \$infile,
	"output=s"   => \$output,
	"dot=s"      => \$outfile,
	"groups=s"   => \$groupfile,
	"stats=s"    => \$statfile,
	"tool_path=s"=> \$tool_path
);


die $usage
  if ( !$infile);


my $HAPLOPHYLE_EXE = "java -Xmx2048m -jar NetworkCreator_fat.jar";
if ($tool_path){
	$HAPLOPHYLE_EXE = "java -Xmx2048m -jar $tool_path/NetworkCreator_fat.jar";
}
  
	
my $out_png = "network.png";

my $command = "$HAPLOPHYLE_EXE -in $infile -out $outfile";
system($command);



my %groups;
my %groups2;
my %hash;
my %haplosize;
if ($groupfile && $statfile){
	open(G,$groupfile);
	while(<G>){
		my $line = $_;$line=~s/\n//g;$line=~s/\r//g;
		my ($ind,$group) = split(";",$line);
		if ($group =~/\w+/ && $ind=~/\w+/){
			$groups{$group}.=$ind.",";	
			$groups2{$ind} = $group;
		}
	}
	close(G);

	open(S,$statfile);
	while(<S>){
		if (/^(haplo\d+):(\d+):(.*)/){
			my $haplo_num = $1;
			my $nb_haplo = $2;
			my $inds = $3;
			my @inds = split(",",$3);
			foreach my $ind(@inds){
				my ($indname,$rank) = split("_",$ind);
				my $group = $groups2{$indname};
				$hash{$haplo_num}{$group}++;
				$haplosize{$haplo_num}++;
			}
		}
	}
	close(S);	
}

my $nb_groups = scalar keys(%groups);
my @colors = ("#ed292a","#ed292a","#82ABA0","#2255a6","#6ebe43","#e76599","#662e91","#c180ff","#ea8b2f","#fff100","#666666","#01ffff","#bfbfbf","#2ac966","#666666");
my $pie_block = "";
my %correspondence_groups;
my $i = 0;
foreach my $group(keys(%groups)){
	$i++;
	$correspondence_groups{$i} = $group;
	$pie_block .= "'pie-$i-background-color': '$colors[$i]',\n";
	$pie_block .= "'pie-$i-background-size': 'mapData(group$i, 0, 10, 0, 100)',\n";
}
open(JSON_CYTOSCAPE,">$output");
my $json = "{\"elements\": {\"nodes\": [";
my $done = 0;
open(OUTFILE,"$outfile");
while(<OUTFILE>){
                        if (/(^\w+)\s\[.*width=([\d\.]+),/){
                                my $node = $1;
                                my $size = $2;
                               	my $ref_hash = $hash{$node};
                                if ($ref_hash){
                                        my %hash2 = %$ref_hash;
                                        my $s = scalar keys(%hash2);
					$json.= "{ \"data\": { \"id\": \"$node\", \"width\": $size";
                                        for (my $i = 1; $i <= $nb_groups; $i++){
						my $group = $correspondence_groups{$i};
						my $ratio = 0;
						if ($haplosize{$node} > 0 && $hash{$node}{$group} > 0){	
	                                                $ratio = ($hash{$node}{$group}/$haplosize{$node}) * 10;
						}
						$json .= ", group$i: $ratio";
                                        }
					$json.= " } },\n";
                                }
                                else{
					$json.= "{ \"data\": { \"id\": \"$node\", \"width\": $size} },\n";
                                }
                        }
                        if (/(\w+) -- (\w+)/){
                                if ($done == 0){
                                        $done = 1;
					chop($json);
					chop($json);
					$json .= "],\n";
					$json .= "\"edges\": [\n";
                                }
                                $done = 1;
				$json.= "{ \"data\": { \"id\": \"$1$2\", \"weight\": 1, \"source\": \"$1\", \"target\": \"$2\"} },\n";
                        }
}
chop($json);
chop($json);
$json.="]}}";
print JSON_CYTOSCAPE $json;
close(JSON_CYTOSCAPE);

