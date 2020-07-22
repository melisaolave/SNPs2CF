# This tutorial shows how to convert a phased diploid vcf file into phylip format.
# NOTE: The function is compatible with the output vcf from SHAPEIT (https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html)

# you will need:
  # a phased vcf file (only diploid allowed). 
  # the library pegas installed. 
    # install pegas (if you don't have it already)
    install.packages("pegas", repos="http://R-Forge.R-project.org");

### Tutorial starts here
# load library
library(pegas)

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

# To convert it into phylip format use the function phased.vcf2phylip:
phased.vcf2phylip(vcf.name="myvcf.vcf", total.SNPs=90, output.name=NULL)

# a new file has been created with the name myvcf.vcf.phy. Check your example folder.
  