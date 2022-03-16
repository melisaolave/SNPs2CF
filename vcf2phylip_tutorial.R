# This tutorial shows how to convert a phased diploid vcf file into phylip format.
# NOTE: The function is compatible with the output vcf from SHAPEIT (https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html).

# you will need:
  # a phased vcf file (only diploid allowed). 
  # the library pegas installed. 
    # install pegas (if you don't have it already)
    install.packages("pegas", repos="http://R-Forge.R-project.org");

### Tutorial starts here
# load libraries
library(pegas)
library(foreach)
library(doMC)

# Then, copy the folder in www.github.com/melisaolave/SNPs2CF to your working directory
# load the functions using source(). This will load the SNPs2CF(), as well as other required internal functions.
# replace MYPATH for your path to the SNPs2CF folder
source("MYPATH/SNPs2CF/functions.R");

# set your working directory (replace MYPATH for your path to the SNPs2CF folder)
setwd("/MYPATH/SNPs2CF/examples/");

# the example folder holds the file myvcf.vcf as an example of a phased vcf file (461 individuals and 90 SNPs)
# you can check it out using:
vcf <- read.vcf("myvcf.vcf", to=90)
vcf

# To convert it into phylip format use the function vcf2phylip:
vcf2phylip(vcf.name="myvcf.vcf", total.SNPs=90, output.name=NULL);

# a new file has been created with the name myvcf.vcf.phy. Check your example folder.

# To randomly assign a base one allele (for non-phased a vcf) use random.phase=T
# NOTE: SNPs2CF assumes SNPs are unlinked. If your dataset does not violate this assumption (e.g. SNPs separated along the chromosome so that linkage is not an issue), then it is only OK to use random.phase = T. For RADseq or target sequencing, usually sampling one SNP per locus would fit this assumption well. For whole genome sequences, you could sampled one SNP every X Kb (e.g. 10Kb, but this depends on recombination map of your specific taxonomic group).
vcf2phylip(vcf.name="myvcf.vcf", total.SNPs=90, random.phase=T);

# To run on parallel use cores
vcf2phylip(vcf.name="myvcf.vcf", total.SNPs=90, cores=8);

# To replace missing data coded as "." in the vcf by "?", use replace.missing = T (as default)
vcf2phylip(vcf.name="myvcf.vcf", total.SNPs=90, replace.missing=T);
  
