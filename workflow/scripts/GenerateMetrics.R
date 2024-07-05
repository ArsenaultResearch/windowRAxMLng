### Clear the Workspace
rm(list=ls()) 

### Load in Necessary libraries
library(ape)
library(SNPRelate)
library(gdsfmt)
library(SeqArray)
library(readr)
library(dplyr)
library(tidyr)

## Interpret Inputs
args <- commandArgs(trailingOnly = TRUE)
window_list <- read_csv(args[1]) ## Load in the window list 
pheno_labels <- read_delim(args[2]) ## load in the phenotype labels as a csv

## Troubleshooting
# setwd("~/Desktop/Arsenault_Research/Bioinformatics_Pipelines/windowRAxMLng/")
# window_list <- read.csv("~/Desktop/Arsenault_Research/Bioinformatics_Pipelines/windowRAxMLng/config/windows_test.csv")
# pheno_labels <- read_delim("~/Desktop/Arsenault_Research/Bioinformatics_Pipelines/windowRAxMLng/config/samples.csv")

pheno_labels_mat <- as.matrix(pheno_labels[,2])
dimnames(pheno_labels_mat) <- list(pheno_labels$Sample_ID,"Caste")
ergatoid.list <- pheno_labels[pheno_labels$Caste=="ergatoid",1]
ergatoid.list <- ergatoid.list$Sample_ID
gynomorph.list <- pheno_labels[pheno_labels$Caste=="gynomorph",1]
gynomorph.list <- gynomorph.list$Sample_ID

## Load analysis functions to run within loop
get_distance <- function(input_tree,ID){
  DistMat <- cophenetic(input_tree)
  
  e2e_dist <- DistMat[row.names(DistMat) %in% ergatoid.list,colnames(DistMat) %in% ergatoid.list]
  e2e_dist[!lower.tri(e2e_dist)] <- NA # remove diagonal and redundant values
  e2e_dist <- data.frame(e2e_dist)
  e2e_long <- pivot_longer(e2e_dist,values_drop_na = T,cols=colnames(e2e_dist))
  
  e2g_dist <- DistMat[row.names(DistMat) %in% gynomorph.list,colnames(DistMat) %in% ergatoid.list]
  e2g_long <- pivot_longer(data.frame(e2g_dist),values_drop_na = T,cols=colnames(data.frame(e2g_dist)))
  
  t.out <- t.test(e2g_long$value,e2e_long$value,paired = F)
  
  mean_frac <- mean(e2e_long$value)/mean(e2g_long$value)
  
  result <- data.frame("tree.mean.ratio"=mean_frac,"tree.t.stat"=t.out$statistic, "tree.t.pval"=t.out$p.value)
  
  return(result)
}

one_snip_test <- function(input_tree){
  l = subtrees(input_tree)
  subtree_res <- data.frame("subtree"= 1:length(l),"metric"=rep(NA,length(l)))
  for (i in 1:length(l)){
    subtree_res[i,2] <- (sum(l[[i]]$tip.label %in% ergatoid.list) / length(ergatoid.list)) - (sum(l[[i]]$tip.label %in% gynomorph.list) / length(gynomorph.list))
  }
  return(max(subtree_res[,2]))
}

make_marked_tree <- function(input_tree,res_folder,ID) {
  pdf(file=paste(res_folder,ID,".markedTree.pdf",sep=""),width = 8,height = 6,onefile = T,useDingbats = F)
  plot(input_tree, x.lim = 1.5,cex=0.5,type="tidy")
  phydataplot(pheno_labels_mat, tree, "m", border = "white", offset = 0.15, width = 0.05)
  dev.off()
}

vcf_dist <- function(ID,res_folder){
  vcf.fn <- paste("./results/",ID,".vcf",sep="")
  gds.fn <- paste("./results/",ID,".gds",sep="")
  
  snpgdsVCF2GDS(vcf.fn, gds.fn, method="biallelic.only") ## vcf.fn is the filename of the vcf
  genofile <- snpgdsOpen(gds.fn)
  pca <- snpgdsPCA(genofile, num.thread=2,autosome.only=FALSE)
  showfile.gds(closeall=TRUE)
  
  pc.percent <- pca$varprop*100
  pc.percent <- round(pc.percent,2)
  PCgt10 <- sum(pc.percent>=10,na.rm = T)
  
  tab <- data.frame(sample.id = pca$sample.id,
                    EV1 = pca$eigenvect[,1],    # the first eigenvector
                    EV2 = pca$eigenvect[,2],    # the second eigenvector
                    stringsAsFactors = FALSE)
  tab$pop <- NA
  tab[tab$sample.id %in% ergatoid.list,c("pop")] <- 1
  tab[tab$sample.id %in% gynomorph.list,c("pop")] <- 2
  
  pdf(file=paste(res_folder,ID,".pca.pdf",sep=""))
  plot(tab$EV1, tab$EV2, col=as.integer(tab$pop), xlab=paste("EV1 - ",pc.percent[1],sep=""), ylab=paste("EV2 - ",pc.percent[2],sep=""))
  dev.off()
  
  tab_mat <- as.matrix(tab[,c(2,3)])
  row.names(tab_mat) <- tab$sample.id
  tab_mat_dist <- as.matrix(dist(tab_mat))
  
  e2e_dist <- tab_mat_dist[row.names(tab_mat_dist) %in% ergatoid.list,colnames(tab_mat_dist) %in% ergatoid.list]
  e2e_dist[!lower.tri(e2e_dist)] <- NA # remove diagonal and redundant values
  e2e_dist <- data.frame(e2e_dist)
  e2e_long <- pivot_longer(e2e_dist,values_drop_na = T,cols=colnames(e2e_dist))
  
  e2g_dist <- tab_mat_dist[row.names(tab_mat_dist) %in% gynomorph.list,colnames(tab_mat_dist) %in% ergatoid.list]
  e2g_long <- pivot_longer(data.frame(e2g_dist),values_drop_na = T,cols=colnames(data.frame(e2g_dist)))
  
  t.out <- t.test(e2g_long$value,e2e_long$value,paired = F)
  
  mean_frac <- mean(e2e_long$value)/mean(e2g_long$value)
  
  result <- data.frame("pca.mean.ratio"=mean_frac,"pca.t.stat"=t.out$statistic, "pca.t.pval"=t.out$p.value,"PCgt10"=PCgt10)

  return(result)
}

## Overarching loop to go through all windows and generate results data.frame

Output <- window_list
namevector <- c("tree.mean.ratio","tree.t.stat","tree.t.pval","one.snip.metric","pca.mean.ratio","pca.t.stat","pca.t.pval","pcaPCgt10")
Output[,namevector] <- NA
for (i in 1:nrow(window_list)){
  tree.name <- paste("./results/",window_list[i,4],".raxml.bestTree",sep="")

  tree <- read.tree(tree.name)
  
  ## Compute Metrics
  dist_result <- get_distance(tree,window_list[i,4])
  Output[i,5:7] <- dist_result
  onesnip_result <- one_snip_test(tree)
  Output[i,8] <- onesnip_result
  vcf_dist_result <- vcf_dist(window_list[i,4],"./results/")
  Output[i,9:12] <- vcf_dist_result
  make_marked_tree(tree,"./results/",window_list[i,4])
}

write_delim(x = Output,file = "./results/metrics.tsv",delim="\t",col_names = T)



