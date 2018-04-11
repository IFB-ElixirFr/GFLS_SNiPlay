#!/usr/bin/perl

use strict;
use Getopt::Long;
use Bio::SeqIO;


my $usage = qq~Usage:$0 <args> [<opts>]
where <args> are:
    -g, --geno         <Genotype input>
    -i, --info         <SNP information. Genome position.>
    -p, --pheno        <Phenotype input>
    -o, --out          <output name>
    -d, --directory    <directory for MLMM R libraries>
    -s, --step_number  <number of steps. Maximum: 20. Default: 10>
    -m, --method       <Method: mbonf or extBIC. Default: mbonf>
~;
$usage .= "\n";

my ($geno,$map,$pheno,$out,$dir,$steps,$method);


GetOptions(
	"geno=s"      => \$geno,
	"info=s"      => \$map,
	"pheno=s"     => \$pheno,
	"out=s"       => \$out,
	"dir=s"       => \$dir,
	"steps=s"     => \$steps,
	"method=s"    => \$method
);


die $usage
  if ( !$geno || !$map || !$pheno || !$out || !$dir || !$steps || !$method);

my $max_steps = 10;
my $plot_opt = "mbonf";
if ($method && $method ne 'mbonf' && $method ne 'extBIC')
{
	print "Aborted: Method must be mbonf or extBIC.\n";
	exit;
}
else
{
	$plot_opt = $method;
}
if ($steps && $steps !~/\d+/ && $steps > 20 && $steps < 2)
{
	print "Aborted: Number of steps must be greater than 2 and lower than 20.\n";
	exit;
}
else
{
	$max_steps = $steps;
}
	

my $chunk = 2;

my $RSCRIPT_EXE = "Rscript";
my $R_DIR = $dir;


my $head_trait = `head -1 $pheno`;
my @headers_traits = split(/\t/,$head_trait);
my $trait_name = $headers_traits[1];


open( my $RCMD, ">rscript" ) or throw Error::Simple($!);



	
print $RCMD "Y_file <- \"" . $pheno . "\"\n";
print $RCMD "X_file <- \"" . $geno . "\"\n";
if($map)
{
	print $RCMD "map_file <- \"$map\"\n";
	print $RCMD "map <- read.table(map_file, sep = \"\\t\", header = T)\n";
}
print $RCMD "mlmm_data = list()\n";
print $RCMD "mlmm_data\$chunk <- $chunk\n";
print $RCMD "mlmm_data\$maxsteps <- $max_steps\n";
print $RCMD "genot <- read.table(X_file, sep = \"\\t\", header = T)\n";
print $RCMD "genot_mat <- as.matrix(genot[, 2:ncol(genot)])\n";
print $RCMD "rownames(genot_mat) <- genot\$Ind_id\n";
	
print $RCMD "phenot <- read.table(Y_file, sep = \"\\t\", header = T)\n";
	
	
	
# missing data imputation
print $RCMD "genot_imp <- genot_mat\n";
print $RCMD "average <- colMeans(genot_imp, na.rm = T)\n";
print $RCMD "for (i in 1:ncol(genot_imp)){genot_imp[is.na(genot_imp[,i]), i] <- average[i]}\n";
	
# kinship matrix computation
print $RCMD "average <- colMeans(genot_imp, na.rm = T)\n";
print $RCMD "stdev <- apply(genot_imp, 2, sd)\n";
print $RCMD "genot_stand <- sweep(sweep(genot_imp, 2, average, \"-\"), 2, stdev, \"/\")\n";
print $RCMD "K_mat <- (genot_stand %*% t(genot_stand)) / ncol(genot_stand)\n";
print $RCMD "write.table(K_mat, '$out.kinship', sep='\\t', dec='.', quote=F, col.names=T, row.names=T)\n";

print $RCMD "source(\"" . $R_DIR. "/mlmm.r\")\n";
print $RCMD "source(\"" . $R_DIR. "/emma.r\")\n";

# mlmm
print $RCMD "mygwas <- mlmm(Y = phenot\$$trait_name, X = genot_imp, K = K_mat, nbchunks=mlmm_data\$chunk, maxsteps=mlmm_data\$maxsteps)\n";
	
# plots
print $RCMD "pdf('$out.pdf')\n";
print $RCMD "plot_step_table(mygwas, \"h2\")\n";
print $RCMD "plot_step_table(mygwas, \"extBIC\")\n";
print $RCMD "plot_step_table(mygwas, \"maxpval\")\n";
print $RCMD "plot_step_RSS(mygwas)\n";
# for (my $j = 1; $j <= ($max_steps - 1); $j++)
# {
	# print $RCMD "plot_fwd_GWAS(mygwas, step = $j, snp_info = map, pval_filt = 0.1)\n";
# }
print $RCMD "plot_opt_GWAS(mygwas, opt = \"extBIC\", snp_info = map, pval_filt = 0.1)\n";
print $RCMD "plot_opt_GWAS(mygwas, opt = \"mbonf\", snp_info = map, pval_filt = 0.1)\n";
#print $RCMD "qqplot_fwd_GWAS(mygwas, nsteps = mlmm_data\$maxsteps)\n";
print $RCMD "qqplot_opt_GWAS(mygwas, opt = \"extBIC\")\n";
print $RCMD "qqplot_opt_GWAS(mygwas, opt = \"mbonf\")\n";

# outputs
print $RCMD "write.table(mygwas\$RSSout, '$out.rss', sep='\\t', dec='.', quote=F, col.names=T, row.names=F)\n";
print $RCMD "write.table(mygwas\$step_table, '$out.steptable', sep='\\t', dec='.', quote=F, col.names=T, row.names=F)\n";

$plot_opt = "\$opt_" . $plot_opt;
print $RCMD "pval = mygwas" . $plot_opt . "\$out\n";
print $RCMD "colnames(pval) = c(\"Marker_name\", \"Pvalue\")\n";
print $RCMD "info_tmp = map\n";
print $RCMD "colnames(info_tmp) = c(\"Marker_name\", \"Chr\", \"Pos\")\n";
print $RCMD "res_asso = pval\n";
print $RCMD qq~
if(exists("info_tmp")){
	res_asso = merge(info_tmp, res_asso, by="Marker_name")
	if( !is.element("Trait", colnames(info_tmp)) ){
		m = matrix(data="traitname", ncol=1, nrow=nrow(res_asso), dimnames=list(c(), c("Trait")))
		res_asso = cbind(m, res_asso)
	}
}
~;
print $RCMD "res_asso = res_asso[order(res_asso[, \"Trait\"], res_asso[, \"Chr\"], res_asso[, \"Pos\"]), ]\n";
print $RCMD "write.table(res_asso, '$out.res_asso', sep='\t', dec='.', quote=F, col.names=T, row.names=F)\n";
close($RCMD);

system("$RSCRIPT_EXE --vanilla rscript");




