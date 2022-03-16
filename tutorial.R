############################################################# SNPs2CF() ###########################################################################
# by Melisa Olave
# Please cite github.com/melisaolave/SNPs2CF (please check for updates in citation)
# For questions or to report issues, refer to molave@mendoza-conicet.gob.ar

############################################################# DESCRIPTION #########################################################################
# This tutorial allows to obtain the concordance factors (CF) calculations from a SNP matrix in phylip format using the SNPs2CF() function.
# Make sure of reading the documentation before.
# When using your own dataset, you can convert your vcf into phylip format using the vcf2phylip (included). See vcf2phylip tutorial form more information
# Version 1.5
#################################### Getting started: Packages installation and loading the functions #############################################
# Make sure to install the doMC and foreach libraries. If they are not yet installed, run:
install.packages("foreach", repos="http://R-Forge.R-project.org");
install.packages("doMC", repos="http://R-Forge.R-project.org");

# Then, copy the folder in www.github.com/melisaolave/SNPs2CF to your working directory
# load the functions using source(). This will load the SNPs2CF(), as well as other required internal functions.
  # replace MYPATH for your path to the SNPs2CF folder
source("MYPATH/SNPs2CF/functions_v1.5.R");

#################################### Using SNF2CF() - 1 individual per species ###################################################################
# set your working directory (replace MYPATH for your path to the SNPs2CF folder)
setwd("/MYPATH/SNPs2CF/examples/");

# In the examples folder there is a SNP phylip matrix from simulated data, named as 5taxa-30K_SNPs.phy.
# You can take a look in a text editor such as WordPad or TextWrangler. The matrix has 5 species and 30,000 SNPs

# Here, we run SNPs2CF. We will set boostrap=FALSE and max.SNPs = 1000 just to use the first 1,000 SNPs, 
  # so the calculations will not take that long. If your computer is still too slow, you can reduce the number of max.SNPs to make things go faster. 
  # NOTE: when using your real data, it is good to take advantage of the largest number of SNPs possible. So, unless there is a clear
  # reason for it, DO NOT use max.SNPs (either don't write this argument or set max.SNPs = NULL)
output <- SNPs2CF(seqMatrix="5taxa-30K_SNPs.phy", max.SNPs = 1000, bootstrap=FALSE, outputName="SNPs2CF_5taxa-1K_SNPs.csv", save.progress=FALSE);

# it is possible to take a look to the output
head(output)

# Check your working directory and should find the new file: SNPs2CF_5taxa-2K_SNPs.csv
# For interpretation of the output check the documentation (Output section)

# Let's now try to include bootstrap and compare the output with the one generated before
output <- SNPs2CF(seqMatrix="5taxa-30K_SNPs.phy", max.SNPs = 1000, bootstrap=TRUE, outputName="SNPs2CF_5taxa-1K_SNPs-bootstrap.csv", save.progress=FALSE);

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
output <- SNPs2CF(seqMatrix="5taxa-2ind-30K_SNPs.phy", ImapName="Imap.txt", between.sp.only = TRUE, max.SNPs = 100, bootstrap=FALSE,
                  outputName="SNPs2CF_5taxa-2ind-1K_SNPs.csv", save.progress=FALSE)

# Check your working directory and should find the new file: SNPs2CF_5taxa-2ind-1K_SNPs.csv

# Let's now try out the between.sp.only = FALSE (default) 
  # IMPORTANT: the number of quartets rises now to 2,100. So, here we reduce the max.SNPs to only 50.
output <- SNPs2CF(seqMatrix="5taxa-2ind-30K_SNPs.phy", ImapName="Imap.txt", between.sp.only = FALSE, max.SNPs = 50, 
                   bootstrap=FALSE, outputName="5taxa-2ind-1K_SNPs.phy-allQuart.csv", save.progress=FALSE);

# Check your working directory and should find the new file: 5taxa-2ind-1K_SNPs.phy-allQuart.csv

#################################### Using SNF2CF() - Subsampling quartets when having multiple individuals per species ################################
# set your working directory (replace MYPATH for your path to the SNPs2CF folder)
setwd("/MYPATH/SNPs2CF/examples/");

# When having multiple individuals per species it is possible to subsample the number of quartets.
# Check the documentation for details about how this works

# Use n.quartets to subsample quartets. Here we will only sample 2 individual quartets per species quartet, by setting n.quartets = 2 
output <- SNPs2CF(seqMatrix="5taxa-2ind-30K_SNPs.phy", ImapName="Imap.txt", between.sp.only = TRUE, max.SNPs = 100, bootstrap=FALSE, n.quartets = 2,
                    outputName="SNPs2CF_5taxa-2ind-1K_SNPs-2quart.csv", save.progress=FALSE);

# Check your working directory and should find the new file: SNPs2CF_5taxa-2ind-1K_SNPs-2quart.csv

#################################### Using SNF2CF() - running in multple cores is parallel #############################################################
# set your working directory (replace MYPATH for your path to the SNPs2CF folder)
setwd("/MYPATH/SNPs2CF/examples/");

# It is possible to performe the calculation of each species quartet in a separate core, to improve time consumption.

# Use cores to set the number of cores. Here we will run the calculations in 2 cores
  # WARNING: before running make sure you have 2 cores available. If there are not enough resources, your computer might crash.
output <- SNPs2CF(seqMatrix="5taxa-30K_SNPs.phy", between.sp.only = TRUE, max.SNPs = 1000, bootstrap=FALSE, save.progress=FALSE,
                  outputName="SNPs2CF_5taxa-1K_SNPs.csv", 
                  cores = 2);

# Check your working directory and should find the new file: SNPs2CF_5taxa-1K_SNPs.csv


################################### Saving progress, continuing calculations and combining temporal CF tables ##########################################
# Here, I described some useful utility to apply for cases of unfinished runs (function crashes, the computer goes off, etc).

# set your working directory (replace MYPATH for your path to the SNPs2CF folder)
setwd("/MYPATH/SNPs2CF/examples/");

# If save.progress = T when running SNPs2CF (default), then a temporal folder is created and the calculation for each species quartet is saved separtely
# Try the code below, to save the progress
output <- SNPs2CF(seqMatrix="5taxa-2ind-30K_SNPs.phy", ImapName="Imap.txt", between.sp.only = TRUE, max.SNPs = 100, bootstrap=FALSE, n.quartets = 2,
                  outputName="SNPs2CF_5taxa-2ind-1K_SNPs-2quart.csv", save.progress=TRUE);

# There is a temporal folder that is created on the working directory. It contains a separated table for each species quartet that were created while 
# the loop was making progress.

# Lets suppose that the electricity went off in your computer and only the CF for the first 3 quartets were calculated. Then, it is possible to use 
# the starting.sp.quartet element to let know the SNPs2CF function should start from the quartet number 4.
output <- SNPs2CF(seqMatrix="5taxa-2ind-30K_SNPs.phy", ImapName="Imap.txt", between.sp.only = TRUE, max.SNPs = 100, bootstrap=FALSE, n.quartets = 2,
                  outputName="SNPs2CF_5taxa-2ind-1K_SNPs-2quart.csv", save.progress=TRUE,
                  starting.sp.quartet=4);

# It is then possible to combine all the partial tables into one, using the function combind.CF.table() as follow:
# if folder.names = NULL, then the function will append all the files within all folders matching the pattern "temp"
# if you want to target a specific folder, use a character vector to name them using folder.names. E.g. folder.names=c("temp_11-25-03", "temp_11-11-03);

combined.output <- combine.CF.table(folder.names=NULL, pattern=".temp.csv$", output="SNPs2CF.csv")

# Then the partial runs were combined and saved in a single table on your working directory
# check it out also:
combined.output

################################### Loading the SNPs2CF output into PhyloNeworks #######################################################################
# Once the CF table is obtained, it is possible to continue with the previously available PhyloNetworks pipeline.

# The output table is in the same format than the one required by PhyloNetworks, and can be loaded into julia (see https://julialang.org):
# julia> using PhyloNetworks
# julia> CF = readTableCF("SNPs2CF.csv")

# if you are using between.sp.only = FALSE, then the CF table has more than one species per line and in such case PhyloNetworks requires one additional step to read the table. Just do:

#julia> using PhyloNetworks
#julia> CFtable = CSV.read("SNPs2CF.csv", copycols=true)
#julia> CF = readTableCF!(CFtable)


# For further steps and phylogenetic network reconstruction, continue with the tutorials
# on the website http://crsl4.github.io/PhyloNetworks.jl/latest/man/snaq_plot/#Network-Estimation-1
