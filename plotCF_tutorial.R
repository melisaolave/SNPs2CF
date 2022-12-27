# Tutorial to plot CF bars from for quartet comparisons

#wd = define working directory path, default current directory

#CF.table.name = define concordance factor table name, e.g. SNPs2CF.txt

#target.quartet= only if want to target specific quartets. Character string. Use comma, no space:
                #e.g. target.quartet=c("sp1,sp2") will only plot quartets containing both sp1 and sp2
                #e.g. target.quartet=c("sp1,sp2" , "sp3,sp4") will plot quartets containing sp1 and sp2, as well as all quartets with sp3 and sp4
                #e.g. target.quartet=c("sp1,sp2,sp3,sp4") will only plot quartet for with all sp1, sp2, sp3 and sp4 
                #by default all quartets are plot

#col = define color for bars, default = gray

#plot.stats = TRUE, if want to plot statistical results (t-test) using asterisks

#p.lines.cex = width of lines for p values

#y.line.adj = change to adjust vertical position of stat line comparisons, default 0.02


##### tutorial starts here #####
# set working directory
setwd("/my/path/goes/here")

#load functions
source("/my/path/to/functions/function_v1.6.R")

#run plotCF for all quartet
plotCF(wd=getwd(), CF.table.name="SNPs2CF.csv", col="gray", 
      plot.stats=TRUE, asterisk.cex=0.85, p.lines.cex=1, y.line.adj=0.02)
      
#if want to plot many bar charts in a single plot, for example 2x5 plot, you can use:
rows <- 2 # number of rows for plotting
columns <- 5 #number of columns for plotting

plot.new()
mymatrix <- matrix(1:(rows*columns), nrow=rows, ncol=columns, byrow=T);
layout(mymatrix);
par(oma = c(6, 2, 0, 2)); #outer margins
par(mar = c(4, 2, 2, 2)); #inner margins


plotCF(wd=getwd(), CF.table.name="SNPs2CF.csv", 
       col="gray", 
      plot.stats=TRUE, asterisk.cex=0.85, p.lines.cex=1, y.line.adj=0.02)

# save plot using
dev.copy(pdf, "CFplot.pdf", width=8, height=6);
dev.off();
