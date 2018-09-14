# SNPs2CF
An R function to compute Concordance Factors from SNP datasets

https://github.com/melisaolave/SNPs2CF
by Melisa Olave, University of Konstanz (Germany)

LICENSE
This software is provided "as is" without warranty of any kind. In no event shall the author be held responsible for any damage resulting from the use of this software. The program package, including the R functions, documentation, tutorials and examples, is distributed free of charge.

DESCRIPTION
This is an R function that performs the quartet-level concordance factor (CF) calculations from single nucleotide polymorphism (SNP) matrices. 
The output generated can later be loaded into the PhyloNetworks julia package (Solis-Lemus et al. 2017; Molecular biology and evolution, 34(12), 3292-3298) to estimate phylogenetic networks. 
Information about PhyloNetworks can be found in its author’s website https://github.com/crsl4/PhyloNetworks.jl

CITATION (Please check for citation updates)
Our pipeline was presented at the II Joint Congress on Evolutionary Biology 2018 (Montpellier, France). Note that it is not yet published on a formal scientific journal but, because of the high demand received, we decided to make this function available before publication.
Please cite:
Olave M. & Meyer A. Tree thinking vs network thinking: a new approach to reconstruct phylogenetic networks from SNP datasets applied to study the rapidly speciating crater lake cichlids from Nicaragua. II Joint Congress on Evolutionary Biology 2018 (Montpellier, France). Symposium: New approaches to Phylogenomics.

REQUIREMENTS
1. R program: https://www.r-project.org
2. R packages: foreach and doMC.

To report issues refer to: melisa.olave@uni-konstanz.de
