install_github("ThomasKAtkins/SNPediaR")
library(SNPediaR)
library(dplyr)
library(readr)

getSnpediaPages <- function(snps) {
  # a wrapper for the SNPediaR getPages function
  # returns NULL in case of error
  out <- tryCatch(
    getPages(snps),
    error = function(e) {
      return(NULL)
    }
  )
  return (out)
}

if (file.exists("snp_df.csv")) {
  print("Warning: snp_df.csv already exists. Please delete before continuing.")
}
if (file.exists("geno_df.csv")) {
  print("Warning: geno_df.csv already exists. Please delete before continuing.")
}

# file from 23andMe
args = commandArgs(trailingOnly=TRUE)
all_snps <- read.table(args[1], sep="\t", header=T)

n <- 1000 #how many SNPs to scrape at a time
is <- seq(1, nrow(all_snps), n)
text <- c()

# set up an empty data frame for genotypes
geno_df <- data.frame(gt1=character(),gt2=character(),gt3=character())

# scrape pages in batches of n
for(i in is) {
  print(i/nrow(all_snps)) # print progress
  pages <- NULL
  # when the call fails, it will return NULL, so we retry until it works
  while(is.null(pages)) {
    pages <- getSnpediaPages(all_snps$rsid[i:(i+n-1)])
  }
  
  # subset all non-NULL pages (of the n we're trying to get)
  pg <- pages[pages!="NULL"]
  # keep going if any pages are non-NULL
  if (length(pg) > 0) {

    df <- data.frame(t(sapply(pg, extractSnpTags))) #from SNPediaR
    df <- df[!is.na(df$rsid),] #only take pages with rsids that aren't NA

    if (nrow(df)==0) {
      next #skip to next batch
    }

    # get the names of the genotype pages to scrape
    df$pg1 <- paste("Rs",df$rsid, df$geno1, sep="")
    df$pg2 <- paste("Rs",df$rsid, df$geno2, sep="")
    df$pg3 <- paste("Rs",df$rsid, df$geno3, sep="")
    
    # scrape the genotype pages
    geno1s <- NULL
    while(is.null(geno1s)) {
      geno1s <- getSnpediaPages(df$pg1)
    }
    geno2s <- NULL
    while(is.null(geno2s)) {
      geno2s <- getSnpediaPages(df$pg2)
    }
    geno3s <- NULL
    while(is.null(geno3s)) {
      geno3s <- getSnpediaPages(df$pg3)
    }
    
    # replace NULLs with the string "none"
    geno1s[geno1s=="NULL"] <- "none"
    geno2s[geno2s=="NULL"] <- "none"
    geno3s[geno3s=="NULL"] <- "none"
    
    # replace newlines with '$' for easier parsing
    g1s <- sapply(geno1s, (function (x) gsub("\\n","$", x)))
    g2s <- sapply(geno2s, (function (x) gsub("\\n","$", x)))
    g3s <- sapply(geno3s, (function (x) gsub("\\n","$", x)))

    # convert to genotype data frame
    g1s <- as.data.frame(g1s)
    g2s <- as.data.frame(g2s)
    g3s <- as.data.frame(g3s)

    # rename the rows of the data frame
    rownames(g1s) <- sapply(rownames(g1s), (function (x) sub("(Rs\\d+)\\(.+\\)","\\1", x)))
    rownames(g2s) <- sapply(rownames(g2s), (function (x) sub("(Rs\\d+)\\(.+\\)","\\1", x)))
    rownames(g3s) <- sapply(rownames(g3s), (function (x) sub("(Rs\\d+)\\(.+\\)","\\1", x)))
    
    # sort the data frames by the renamed rows so that they all sort in the same
    # order
    g1s <- g1s[order(rownames(g1s)),]
    g2s <- g2s[order(rownames(g2s)),]
    g3s <- g3s[order(rownames(g3s)),]

    # create a combined data frame with the genotypes for each SNP
    geno_df <- data.frame(cbind.data.frame(g1s, g2s, g3s))

    # if there are any SNPs with genotype pages
    if (nrow(geno_df)>0) {
      colnames(geno_df) <- c("gt1", "gt2", "gt3")

      # make sure that one genotype has information available
      good <- (geno_df$gt1!="none") + (geno_df$gt2!="none") + (geno_df$gt3!="none")
      geno_df <- geno_df[good>0,]
    }
  }

  # replace all null pages with the empty string
  pages[pages=="NULL"] <- ""
  # replace newlines with '$' for easier parsing
  pages <- sapply(pages, (function (x) gsub("\\n","$", x)))
  # remove quotes so that the csv parser doens't get angry
  pages <- gsub("\"","'", pages)
  geno_df$gt1 <- gsub("\"","'", geno_df$gt1)
  geno_df$gt2 <- gsub("\"","'", geno_df$gt2)
  geno_df$gt3 <- gsub("\"","'", geno_df$gt3)
  # create a df of SNPs parsed and remove all without any text
  snp_df <- data.frame(page=pages)
  snp_df <- snp_df[snp_df$page!="",]
  colnames(snp_df) <- c("text")
  # append to the SNP data frame
  write.table(snp_df, "snp_df.csv", append=T, sep=",", col.names = !file.exists("snp_df.csv"), row.names=F)
  # append to the genotype data frame
  write.table(geno_df, "geno_df.csv", append=T, sep=",", col.names = !file.exists("geno_df.csv"), row.names=F)
}
