SNPs2CF: An R function to compute Concordance Factors from SNP datasets
https://github.com/melisaolave/SNPs2CF
by Melisa Olave, University of Konstanz (Germany)

LICENSE
This software is provided "as is" without warranty of any kind. In no event shall the author be held responsible for any damage resulting from the use of this software. The program package, including the R functions, documentation, tutorials and examples, is distributed free of charge.

DESCRIPTION
This is an R function that performs the quartet-level concordance factor (CF) calculations from single nucleotide polymorphism (SNP) matrices. 
The output generated can later be loaded into the PhyloNetworks julia package (Solis-Lemus et al. 2017; Molecular biology and evolution, 34(12), 3292-3298) to estimate phylogenetic networks. 
Information about PhyloNetworks can be found in its author’s website https://github.com/crsl4/PhyloNetworks.jl

CITATION
Olave M. & Meyer A (2020). Implementing Large Genomic SNP Datasets in Phylogenetic Network Reconstructions: A Case Study of Particularly Rapid Radiations of Cichlid Fish. Systematic Biology (in press).

REQUIREMENTS
1. R program: https://www.r-project.org
2. R packages: foreach and doMC.

To report issues refer to: melisa.olave@uni-konstanz.de

CHANGES:
version 1.41:
- A note about reading CF with multiple alleles was added to the documentation and tutorial

version 1.4:
- the phased.vcf2phylip function was added to functions.R, as well as the tutorial vcf2phylip_tutorial.R

version 1.3:
- a warning is saved in the log file if no informative SNPs are found to inform a species quartet when n.quartets = "all"
- bug fix: bootstrap confidence interval returned an error if there were no SNPs informing CFs in one or more pseudoreplicates. Common error when there are only very few informative SNP in the matrix.

version 1.2: 
- the object save.progress was added to the SNPs2CF function, to allow saving progress
- the new function combine.CF.table was added to the functions.R
