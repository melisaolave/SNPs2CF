library(foreach);
library(doMC);

SNPs2CF <- function(wd=getwd(), seqMatrix, 
                    ImapName=NULL, rm.outgroup=F, outgroupSp="outgroup",
                    indels.as.fifth.state=F,  
                    bootstrap=T, boots.rep=100, 
                    outputName="SNPs2CF.csv",
                    n.quartets="all", between.sp.only=F, max.SNPs=NULL, 
                    cores=1){
  registerDoMC(cores);
  setwd(wd);
  options(scipen=999);
  # getting inputs
  matrix <- scan(seqMatrix, what="character", sep="\n");
  if(is.null(ImapName)){
    cat("Not Imap provided, each allele is considered as a different species\n");
    indName <- gsub(pattern="^(.*) +.*$", replacement="\\1", matrix[2:length(matrix)]);
    indName <- gsub(pattern=" ", replacement="", indName);
    Imap <- data.frame(indName, indName);
    colnames(Imap) <- c("traits", "species");
    if(between.sp.only == F){
      cat("Since an Imap was not provided, then between.sp.only = F is ignored\n");
      between.sp.only <- T;
    }
  }else{
    cat("Reading Imap file\n");
    cat("Two columns are expected, first column are the individuals (or alleles) and second column are the species\n");
    Imap <- read.table(ImapName, sep="\t", header=T);
    colnames(Imap) <- c("traits", "species");
  }
  # Checking Imap errors
  if(length(Imap$traits) == 0 | length(Imap$species) == 0){
    stop("Error in Imap file. Make sure the column names are traits and species, tab-separated\n");
  }else if(length(Imap$traits) != length(unique(Imap$traits))){
    stop("Error in Imap file. Each individual has to have a unique name. Check names in trait column\n");
  }
  
  # removing outgroup
  if(rm.outgroup){
    cat("Removing outgroup", outgroupSp, "\n");
    rm.rows <- grep(pattern=paste(outgroupSp, collapse="|"),Imap$species);
    ind.outNames <- as.character(Imap$species[rm.rows]);
    keep.rows <- setdiff(1:nrow(Imap), rm.rows)
    Imap <- Imap[keep.rows,];
    rm.rows <- grep(pattern=paste(ind.outNames, collapse="|"), matrix);
    keep.rows <- setdiff(1:length(matrix), rm.rows);
    matrix <- matrix[keep.rows];
  }
  
  # Matrix into table
  matrix <- matrix[2:length(matrix)]; #removing phylip header
  indMatrix <- gsub(pattern="^(.*) +.*$", replacement="\\1", matrix);
  indMatrix <- gsub(pattern=" +", replacement="", indMatrix)
  seq <- gsub(pattern="^.* +(.*)$", replacement="\\1", matrix)
  table <- data.frame();
  table <- do.call(rbind, strsplit(seq, ""));
  matrixTable <- cbind(indMatrix, table);
  
  # Checking errors in matrix
  if(length(indMatrix) != length(unique(indMatrix))){
    stop("Error in SNP matrix. Each individual has to have a unique name.\n");
  }
  if(length(setdiff(indMatrix, Imap$traits)) > 0){
    cat("WARNING: Some individual names in the SNP matrix were not found in Imap. Only those listed in the Imap file will be used for CF calculations. Individuals not found:\n", paste(setdiff(indMatrix, Imap$traits), " "), "\n");
  }else if(length(setdiff(Imap$traits, indMatrix)) > 0){
    stop("Error: Not all individuals in the SNP matrix were found in the Imap file:", paste(setdiff(Imap$traits, indMatrix), sep=" "), "\n");
  }
  check.indMatrix <- strsplit(indMatrix, ""); # check that all names in individuals make sense (sequence was not matched as part of the name)
  for(ii in 1:length(indMatrix)){
    if(length(check.indMatrix[[ii]]) > 50){ # I guess no one will have a name with more than 50 letters
      stop("Something went wrong with a name in the SNP matrix. Make sure that all sequence names have one or more spaces separating the SNPs. Error reading:\n", indMatrix[ii], "\n")
    }
  }
  
  # getting quartets
  if(between.sp.only){
    spNames <- sort(as.character(unique(Imap$species)));
    n.sp <-  length(spNames);
    cat(n.sp, "species found in Imap:", spNames, "\n");
    allQuartets <- combn(x=1:n.sp, m=4); # getting all possible quartets combinations
  }else{
    if(n.quartets != "all"){
      cat("WARNING: if between.sp.only = F, then n.quartets will be treated as = all\n")
      n.quartets <- "all";
    }
    spNames <- sort(as.character(Imap$species));
    n.sp <-  length(spNames);
    cat(length(unique(spNames)), "individuals found in Imap:", unique(spNames), "\n");
    allQuartets <- combn(x=1:n.sp, m=4); # getting all possible quartets combinations
    if(ncol(allQuartets) > 100000){
      stop("Too many quartets (>100,000). Set between.sp.only = T, and run again\n");
    }
  }
  cat("All possible species quartets for", n.sp, "species:", ncol(allQuartets), "\n");
  
  # checking other arguments
  if(is.null(max.SNPs) == F){
    cat("WARNING: by setting max.SNPs the maximum number of SNPs per quartet will be restricted to: max.SNPs = ", max.SNPs, "\n");
  }
  if(bootstrap){
    cat("bootstrap was set as TRUE, then a total of ", boots.rep, "bootstrap pseudo-replicates will be computed\n")
  }
  if(cores > 1){
    if(n.quartets == "all"){
      cat("n.quartets = all, thus CF will be calculated from all possible quartets\n"); 
    }else{
      cat("A total of ", ncol(allQuartets)*n.quartets, "quartets will be sampled\n"); 
    }
    cat("Getting CF calculations. There is not progress report when running in parallel. This could take several minutes/hours, please be pacient...\n")
  }
  if(indels.as.fifth.state){
    gaps <- "-";
  }else{
    gaps <- NULL;
  }
  
  # starting loop
  indQuartetVec <- NULL;
  quartet.list <- NULL;
  n.resampled <- 0;
  output.table <- data.frame(matrix(ncol=8, nrow=0));
  error.log <- paste(seqMatrix, "-SNPs2CF.err", sep="");
  write(NULL,error.log);
  resampling.err <- FALSE;
  start.time <- Sys.time();
  output.table <- foreach(l=1:ncol(allQuartets), .combine=rbind) %dopar%{
    cat("Working on species quartet:", l, "/", ncol(allQuartets), "\n");
    sampled.sp <- spNames[allQuartets[,l]];
    if(between.sp.only == F){
      if(sampled.sp[1] == sampled.sp[2] & sampled.sp[2] == sampled.sp[3] & sampled.sp[3] == sampled.sp[4]){ # rows corresponding to 4 individuals from the same species are uninformative about between-species relationships, then break loop
        cat("4 individuals from the same species are uninformative about between-species relationships, breaking loop and continuing\n")
        break;
      }
    }
    cat(sampled.sp, "\n");
    if(n.quartets == "all"){
      sp1 <- as.character(Imap[Imap$species == sampled.sp[1],1]);
      sp2 <- as.character(Imap[Imap$species == sampled.sp[2],1]);
      sp3 <- as.character(Imap[Imap$species == sampled.sp[3],1]);
      sp4 <- as.character(Imap[Imap$species == sampled.sp[4],1]);
      all.indQuartet <- get.allQuartet.combn(sp1, sp2, sp3, sp4);
      n.quart <- ncol(all.indQuartet);
    }else{
      n.quart <- n.quartets;
    }
    temp.table <- data.frame(matrix(ncol=8, nrow=0));
    count <- 1;
    n.resampled <- 0;
    ind.sampling <- 0;
    max.boots.split12_34 <- NULL;
    min.boots.split12_34 <- NULL;
    max.boots.split13_24 <- NULL;
    min.boots.split13_24 <- NULL;
    max.boots.split14_23 <- NULL;
    min.boots.split14_23 <- NULL;
    quartet.ok <- T;
    cat("A total of", n.quart, "individual quartets will be sampled within each of the", ncol(allQuartets), "species quartets\n");
    cat("Progress:\n")
    while(count <= n.quart){
      cat("\t", count, "/", n.quart, "of", l, "/", ncol(allQuartets), "\n");
      if(n.quartets != "all"){
        random.indQuartet <- NULL;
        for(i in 1:4){
          random.indQuartet[i] <- sample(x=as.character(Imap[Imap$species == sampled.sp[i],1]), size=1);
        }
        ind.sampling <- ind.sampling+1;
        allQuartets.logic <- F;
      }else{
        allQuartets.logic <- T;
        random.indQuartet <- all.indQuartet[,count];
      }
      if(between.sp.only == F){
        if(length(unique(all.indQuartet[,count])) != 4){
          dont.skip <- F; # if the same individual was sampled twice or more, then do not calculate CF
        }else{
          dont.skip <- T;
        }
      }else{
        dont.skip <- T;
      }
      random.indQuartetVec <- paste(sort(random.indQuartet), collapse="|");
      if((allQuartets.logic | random.indQuartetVec %in% indQuartetVec == F) & dont.skip){ #if the quartet was nos sampled yet
        # CF from new quartet
        rows <- grep(pattern=paste("^",random.indQuartet, "$",sep="",collapse="|"), indMatrix);
        subMatrix <- matrixTable[rows,];
        for(k in 1:4){
          subMatrix[k,1] <- as.character(Imap[Imap$traits == subMatrix[k,1], 2]);
        }
        subMatrix <- subMatrix[order(subMatrix[,1]),];
        split12_34 <- NULL;
        split13_24 <- NULL;
        split14_23 <- NULL;
        for(j in 2:ncol(subMatrix)){
          bp <- subMatrix[,c(1,j)];
          alleles <- unique(bp[,2]);
          if(length(setdiff(alleles, c("A", "C", "T", "G", "a", "c", "t", "g", as.character(0:10), gaps))) == 0  # ignoring sites with ambiguities, missing data and gaps
             & length(alleles) == 2){ # if there are only bases and only biallelics
            split12_34 <- c(split12_34, bp[1,2] == bp[2,2] & bp[3,2] == bp[4,2]);
            split13_24 <- c(split13_24, bp[1,2] == bp[3,2] & bp[2,2] == bp[4,2]);
            split14_23 <- c(split14_23, bp[1,2] == bp[4,2] & bp[2,2] == bp[3,2]);
            if(is.null(max.SNPs) == F){
              if(sum(split12_34, split13_24, split14_23) == max.SNPs){
                break;
              }
            }
          }
        }
        genes <- sum(split12_34, split13_24, split14_23)
        if(is.null(max.SNPs) == F){
          if(genes < max.SNPs){
            cat("max.SNPs not reached, total of genes for this quartet = ", genes, "\n");
          }
        }
        CFVec <- c(sum(split12_34), sum(split13_24), sum(split14_23));
        if(sum(CFVec) != 0){ # only if there was a minimum of one variable site, else, pick other quartet
          indQuartetVec[count] <- random.indQuartetVec;
          temp.table[count,] <- c(sampled.sp, CFVec, genes);
          count <- count+1;
          ### boostrap ####
          if(bootstrap){ # bootstrap only in those quartets with at least one variable site
            boots.split12_34.vec <- NULL;
            boots.split13_24.vec <- NULL;
            boots.split14_23.vec <- NULL;
            cat("\tCalculating CF for", boots.rep, "bootstrap pseudo-replicates... ")
            for(b in 1:boots.rep){
              boots.snps <- sample(1:length(split12_34), length(split12_34), replace=T); # randomly sampling SNPs
              boots.split12_34 <- split12_34[boots.snps];
              boots.split13_24 <- split13_24[boots.snps];
              boots.split14_23 <- split14_23[boots.snps];
              boots.genes <- sum(boots.split12_34, boots.split13_24, boots.split14_23);
              boots.split12_34.vec[b] <- round(sum(boots.split12_34)/boots.genes, 4);
              boots.split13_24.vec[b] <- round(sum(boots.split13_24)/boots.genes, 4);
              boots.split14_23.vec[b] <- round(sum(boots.split14_23)/boots.genes, 4);
            }
            cat("done!\n");
            max.boots.split12_34 <- c(max.boots.split12_34, quantile(boots.split12_34.vec, 0.95));
            min.boots.split12_34 <- c(min.boots.split12_34, quantile(boots.split12_34.vec, 0.05));
            max.boots.split13_24 <- c(max.boots.split13_24, quantile(boots.split13_24.vec, 0.95));
            min.boots.split13_24 <- c(min.boots.split13_24, quantile(boots.split13_24.vec, 0.05));
            max.boots.split14_23 <- c(max.boots.split14_23, quantile(boots.split14_23.vec, 0.95));
            min.boots.split14_23 <- c(min.boots.split14_23, quantile(boots.split14_23.vec, 0.05));
          }
          ################### bootstrap #####
        }else{
          n.resampled <- n.resampled +1;
        }
      }else{
        cat("quartet re-sampled, discarded...\n");
        n.resampled <- n.resampled +1;
      }
      if(n.resampled > 1){
        cat("\tresampling =", n.resampled,"\n");
      }
      if(n.resampled == 100){ # in case the while loop gets stock for ever when randomly sampling
        error <- paste("Warning: Species quartet:", paste(sampled.sp, collapse=" "), "\nwas sampled", ind.sampling, "times only, instead of n.quartets=", n.quartets, "\n");
        prev.error <- scan(error.log, what="character", sep="\n");
        new.error <- c(prev.error, error);
        write(new.error, error.log);
        quartet.ok <- F; # do not try to save this quartet in the table
        break;
      }
    }
      if(quartet.ok){
        sp1 <- rep(as.character(subMatrix[1,1]), length(genes));
        sp2 <- rep(as.character(subMatrix[2,1]), length(genes));
        sp3 <- rep(as.character(subMatrix[3,1]), length(genes));
        sp4 <- rep(as.character(subMatrix[4,1]), length(genes));
        colnames(temp.table) <- c("t1", "t2", "t3", "t4", "CF12_34", "CF13_24", "CF14_23", "genes");
        genes <- c(as.numeric(temp.table[,8]));
        split12_34 <- round(c(as.numeric(temp.table[,5]))/genes, 4);
        split13_24 <- round(c(as.numeric(temp.table[,6]))/genes, 4);
        split14_23 <- round(c(as.numeric(temp.table[,7]))/genes, 4);
        if(bootstrap){
          output.table <- data.frame(sp1, sp2, sp3, sp4, 
                                     split12_34, min.boots.split12_34, max.boots.split12_34,
                                     split13_24, min.boots.split13_24, max.boots.split13_24,
                                     split14_23, min.boots.split14_23, max.boots.split14_23,
                                     genes);
        }else{
          output.table <- data.frame(sp1, sp2, sp3, sp4, split12_34, split13_24, split14_23, genes);
        }
        output.table;
    }
  }
  if(bootstrap){
    colnames(output.table) <- c("t1", "t2", "t3", "t4", 
                                "CF12_34", "CF12_34_lo", "CF12_34_hi",
                                "CF13_24", "CF13_24_lo", "CF13_24_hi",
                                "CF14_23", "CF14_23_lo", "CF14_23_hi",
                                "genes");
  }else{
    colnames(output.table) <- c("t1", "t2", "t3", "t4", "CF12_34", "CF13_24", "CF14_23", "genes");
  }
  write.table(output.table, outputName, sep=",", col.names=T, row.names=F, quote=F);
  prev.error <- scan(error.log, what="character", sep="\n");
  if(length(prev.error) == 0){
    cat("\nWARNING: check the .err file for errors\n");
  }
  end.time <- Sys.time()
  time.taken <- difftime(end.time, start.time, units="secs");
  cat("done!\n\ta total of",nrow(output.table),"quartets were analyzed\n\tTime taken:", time.taken, "seconds\n");
}


get.allQuartet.combn <- function(sp1, sp2, sp3, sp4, stop.max=100000){
  all.indQuartet <- NULL;
  for(m in 1:length(sp1)){
    for(n in 1:length(sp2)){
      for(o in 1:length(sp3)){
        for(p in 1:length(sp4)){
          if(length(unique(c(sp1[m], sp2[n], sp3[o], sp4[p]))) != 4){ # then the same individual is not sampled twice
            break;
          }else{
            indQuartet <- as.matrix(c(sp1[m], sp2[n], sp3[o], sp4[p]));
            all.indQuartet <- cbind(all.indQuartet, indQuartet);
          }
          if(ncol(all.indQuartet) > stop.max){
            stop("Too many individual quartet combinations (>100,000 for a single species quartet). Please select a subset of quartets to explore, changing n.quartet\n");
          }
        }
      }
    }
  }
  cat("sampling all quartets:", ncol(all.indQuartet), "\n");
  return(all.indQuartet);
}
