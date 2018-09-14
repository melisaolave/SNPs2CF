############################################################# SNPs2CF() ###########################################################################
# by Melisa Olave
# Please cite github.com/melisaolave/SNPs2CF (please check for updates in citation)
# For questions or to report issues, refer to melisa.olave@uni-konstanz.de

############################################################# DESCRIPTION #########################################################################
# This tutorial allows to obtain the concordance factors (CF) calculations from a SNP matrix in phylip format using the SNPs2CF() function.
# Make sure of reading the documentation before going with this tutorial.

#################################### Getting started: Packages installation and loading the functions #############################################
# Make sure to install the doMC and foreach libraries. If they are not yet installed, run:
install.packages("foreach", repos="http://R-Forge.R-project.org");
install.packages("doMC", repos="http://R-Forge.R-project.org");

# Then, copy the folder in www.github.com/melisaolave/SNPs2CF to your working directory
# load the functions using source(). This will load the SNPs2CF(), as well as other required internal functions.
  # replace MYPATH for your path to the SNPs2CF folder
source("MYPATH/SNPs2CF/functions.R");

#################################### Using SNF2CF() - 1 individual per species ###################################################################
# set your working directory (replace MYPATH for your path to the SNPs2CF folder)
setwd("/MYPATH/SNPs2CF/examples/");

# In the examples folder there is a SNP phylip matrix from simulated data, named as 5taxa-30K_SNPs.phy.
# You can give it a look in a text editor such as WordPad or TextWrangler. The matrix has 5 species and 30,000 SNPs

# Here, we run SNPs2CF. We will set boostrap=FALSE and max.SNPs = 1000 just to use the first 1,000 SNPs, 
  # so the calculations will not take that long. If your computer is still too slow, you can reduce the number of max.SNPs to make things go faster. 
  # NOTE: when using your real data, it is good to take advantage of the largest number of SNPs possible. So, unless there is a clear
  # reason for it, DO NOT use max.SNPs (either don't write this argument or set max.SNPs = NULL)
SNPs2CF(seqMatrix="5taxa-30K_SNPs.phy", max.SNPs = 1000, bootstrap=FALSE, outputName="SNPs2CF_5taxa-1K_SNPs.csv");

# Check your working directory and should find the new file: SNPs2CF_5taxa-2K_SNPs.csv
# For interpretation of the output check the documentation (Output section)

# Let's now try to include bootstrap and compare the output with the one generated before
SNPs2CF(seqMatrix="5taxa-30K_SNPs.phy", max.SNPs = 1000, bootstrap=TRUE, outputName="SNPs2CF_5taxa-1K_SNPs-bootstrap.csv");

# Check your working directory and should find the new file: SNPs2CF_5taxa-2K_SNPs-bootstrap.csv
# If bootstrap = TRUE then the output also contains two extra columns per each CF for lower and upper limit of the credibility interval.

#################################### Using SNF2CF() - multiple individuals per species #################################################################
# set your working directory (replace MYPATH for your path to the SNPs2CF folder)
setwd("/MYPATH/SNPs2CF/examples/");

# When having multiple individuals per species, an Imap file is required with the individual - species associations.
# Check the Imap.txt file as example and/or go the documentation for details.

# There is an example provided of a SNP matrix with 2 individuals per species named as 5taxa-2ind-30K_SNPs.phy in /SNPs2CF/examples/
  # The number of total quartets increase now from 5 to 80, thus we will set max.SNPs = 100 which is a very small number,
  # but it reduces the time required and, to the end of the tutorial, this is enough to see how it works.
  # Obtain a CF table as follow:
SNPs2CF(seqMatrix="5taxa-2ind-30K_SNPs.phy", ImapName="Imap.txt", between.sp.only = TRUE, max.SNPs = 100, bootstrap=FALSE,
        outputName="SNPs2CF_5taxa-2ind-1K_SNPs.csv")

# Check your working directory and should find the new file: SNPs2CF_5taxa-2ind-1K_SNPs.csv

# Let's now try out the between.sp.only = FALSE (default) 
  # IMPORTANT: the number of quartets rises now to 2,100. So, here we reduce the max.SNPs to only 50.
SNPs2CF(seqMatrix="5taxa-2ind-30K_SNPs.phy", ImapName="Imap.txt", between.sp.only = FALSE, max.SNPs = 50, 
        bootstrap=FALSE, outputName="5taxa-2ind-1K_SNPs.phy-allQuart.csv");

# Check your working directory and should find the new file: 5taxa-2ind-1K_SNPs.phy-allQuart.csv

#################################### Using SNF2CF() - Subsampling quartets when having multiple individuals per species ################################
# set your working directory (replace MYPATH for your path to the SNPs2CF folder)
setwd("/MYPATH/SNPs2CF/examples/");

# When having multiple individuals per species it is possible to subsample the number of quartets.
# Check the documentation for details about how this works

# Use n.quartets to subsample quartets. Here we will only sample 2 individual quartets per species quartet, by setting n.quartets = 2 
SNPs2CF(seqMatrix="5taxa-2ind-30K_SNPs.phy", ImapName="Imap.txt", between.sp.only = TRUE, max.SNPs = 100, bootstrap=FALSE, n.quartets = 2,
        outputName="SNPs2CF_5taxa-2ind-1K_SNPs-2quart.csv");

# Check your working directory and should find the new file: SNPs2CF_5taxa-2ind-1K_SNPs-2quart.csv

#################################### Using SNF2CF() - running in multple cores is parallel #############################################################
# set your working directory (replace MYPATH for your path to the SNPs2CF folder)
setwd("/MYPATH/SNPs2CF/examples/");

# It is possible to performe the calculation of each species quartet in a separate core, to improve time consumption.

# Use cores to set the number of cores. Here we will run the calculations in 2 cores
  # WARNING: before running make sure you have 2 cores available. If there are not enough resources, your computer might crash.
SNPs2CF(seqMatrix="5taxa-30K_SNPs.phy", between.sp.only = TRUE, max.SNPs = 1000, bootstrap=FALSE, outputName="SNPs2CF_5taxa-1K_SNPs.csv", cores = 2);

# Check your working directory and should find the new file: SNPs2CF_5taxa-1K_SNPs.csv

################################### Loading the SNPs2CF output into PhyloNeworks #######################################################################
# Once the CF table is obtained, it is possible to continue with the previously available PhyloNetworks pipeline.

# The output table is in the same format than the one required by PhyloNetworks, and can be loaded into julia (see https://julialang.org):
# julia> using PhyloNetworks
# julia> CF = readTableCF("SNPs2CF.csv")

# For further steps and phylogenetic network reconstruction, continue with the tutorials
# on the website http://crsl4.github.io/PhyloNetworks.jl/latest/man/snaq_plot/#Network-Estimation-1