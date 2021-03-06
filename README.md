
# GFLS (Galaxy For Life Science)  :heavy_minus_sign::sparkles::sparkles::four::herb::microscope::heavy_minus_sign: SNiPlay

This repository concern the **SNiPlay** part of the GFLS project ([website](http://sniplay.southgreen.fr/cgi-bin/home.cgi)).
It is made in collaboration with Alexis Dereeper from the 'Institut de Recherche pour le Développement' (IRD).
All the **XML** files (Wrappers) are made by his team.

#### :file_folder: Content

The repository contain all the wrappers for the workflow SNiPlay.

- admixture [deprecated]
- fastme
- plink
- readseq
- snmf
- snpeff_from_gff_vcf
- tassel5
- vcftools_filter_stats_diversity (vcftoolsFilter,vcftoolsStats,vcftoolsSlidingWindow)
- mlmm

Wrappers for the Haplotype analysis workflow.

- Beagle
- Haplophyle


It contain also wrappers for SNP analysis (in Sniplay repo.).

- AnnotationStatsFromVCF
- check_gwas_inputs
- egglib
- GetHaplotypesFromPhasedVCF
- hapmap2mlmm
- MDSplot
- ped2bed
- PedToFasta
- Rooting
- SNP_density
- VCF2Hapmap


The sniplay3_complete_workflow repository.

~~And the Haplotype_analysis workflow repository.~~


#### :page_with_curl: Publications

Alexis Dereeper, Felix Homa, Gwendoline Andres, Guilhem Sempere, Gautier Sarah, Yann Hueber, Jean-François Dufayard, Manuel Ruiz; SNiPlay3: a web-based application for exploration and large scale analyses of genomic variations, Nucleic Acids Research, Volume 43, Issue W1, 1 July 2015, Pages W295–W300, https://doi.org/10.1093/nar/gkv351

Alexis Dereeper, Stéphane Nicolas, Loïc Le Cunff, Roberto Bacilieri, Agnès Doligez, Jean-Pierre Peros, Manuel Ruiz and Patrice This; SNiPlay: a web-based tool for detection, management and analysis of SNPs. Application to grapevine diversity projects, BMC Bioinformatics, Volume 12, 5 May 2011, Page 134, https://doi.org/10.1186/1471-2105-12-134


####  :pencil: Use tools

These tools are available in a docker container:

- https://quay.io/repository/valentinmarcon/docker-galaxy-sniplay
- https://github.com/ValentinMarcon/docker-galaxy-sniplay
