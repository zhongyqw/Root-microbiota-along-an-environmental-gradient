if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
library(BiocManager)
BiocManager::install("phyloseq")
BiocManager::install("DESeq2")
install.packages("network")
install.packages("sna")
install.packages("tidyverse")
library(devtools)
install_github("umerijaz/microbiomeSeq") 
library(microbiomeSeq)
devtools::install_github("taowenmicro/ggClusterNet")
#update.packages()
library(devtools)
install_github("zdk123/SpiecEasi")
devtools::install_github("zdk123/SpiecEasi")

#--Import the required R package #-------
library(phyloseq)
library(igraph)
library(network)
library(sna)
library(ggplot2)
library(ggClusterNet)
library(tidyverse)

####### Construct and preserve the finseq objects required by ggClusterNet
# Read from file
#RB
metadata = read.table("RB.metadata.txt", header=T, row.names=1, sep="\t", comment.char="", stringsAsFactors = F)
otutab = read.table("RB.core.asv.16S.txt", header=T, row.names=1, sep="\t", comment.char="", stringsAsFactors = F)
taxonomy = read.table("core.tax.df.txt", header=T, row.names=1, sep="\t", stringsAsFactors = F,quote = "")
#otutab<-otutab[,1:84]


# Extract only those ID in common between the two tables
idx = rownames(otutab) %in% rownames(taxonomy)
otutab = otutab[idx,]
taxonomy = taxonomy[rownames(otutab),]

# Import the finseq (PS) object
ps = phyloseq(sample_data(metadata),otu_table(as.matrix(otutab), taxa_are_rows=TRUE), tax_table(as.matrix(taxonomy)))

saveRDS(ps, "ps_RB.rds")

path = "./RB_result-sprcc-type"
dir.create(path)
result = network(ps = ps,N=300, r.threshold=0.6,p.threshold=0.05,
                 label = FALSE,path = path ,group = "type",zipi = TRUE,method = "sparcc")


# Network comparison of all samples
p = result[[1]]
p
# Comparison of network parameters of all samples
data = result[[2]]
plotname1 = paste(path,"/network_all.jpg",sep = "")
ggsave(plotname1, p,width = 22,height = 22)

plotname1 = paste(path,"/network_all.pdf",sep = "")
ggsave(plotname1, p,width = 22,height = 26)

tablename <- paste(path,"/co-occurrence_Grobel_net",".csv",sep = "")
write.csv(data,tablename)