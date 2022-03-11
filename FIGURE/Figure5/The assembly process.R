library(tidyverse)
library(Biostrings)
library(phyloseq)
metadata = read.delim("../data/metadata.txt",row.names = 1)
otutab = read.delim("../data/otutab.txt", row.names=1)
taxonomy = read.delim("../data/Taxonomy.txt", row.names=1)
tree  = read_tree("../data/otus.tree")
rep = readDNAStringSet("../data/otus.fa")

ps = phyloseq(sample_data(metadata),
              otu_table(as.matrix(otutab), taxa_are_rows=TRUE),
              tax_table(as.matrix(taxonomy)), phy_tree(tree),refseq(rep)
)

#- Theme - color etc ------#-------
source("./total_amplicon.R")
#-- Amplicon environment layout

res = theme_my()
mytheme1 = res[[1]]; mytheme2 = res[[2]];  colset1 = res[[3]]; colset2 = res[[4]]; colset3 = res[[5]]; colset4 = res[[6]]


# microbiome data input #------

#-- Sample filtering -- based on information in map files such as grouping and ID #-------
# ps_subsubset_samples(ps0,! ID %in% c("sample1")); ps_sub
# Sequence screening - seven phylum information of rhizosphere or ID information of OTU #---------
ps0 <- ps0 %>%
  # subset_taxa(
  #   Kingdom == id
  # ) %>% 
  phyloseq::filter_taxa(function(x) sum(x ) > 0 , TRUE);ps0


ps = ps0

# # # # # # # # #
library(picante)
library(ape)
library(vegan)
library(FSA)
library(eulerr)
library(grid)
library(gridExtra)
require(minpack.lm)
require(Hmisc)
require(stats4)
library(parallel)
library(ggClusterNet)


psphy = filter_taxa(ps, function(x) sum(x ) > 2000 , TRUE);psphy


res1path = "./"
phypath = paste(res1path,"/Phylogenetic_analyse_spacies/",sep = "")
dir.create(phypath)

map = sample_data(ps)
n = map$Group %>% unique() %>%
  length()
n


#-- Neutral model #-----
source("./phylo_Micro\\neutralModel.R")
result = neutralModel(ps = psphy,group  = "Group",ncol = 4)
#-- Merge charts
p1 =  result[[1]]
p1


FileName <- paste(phypath,"1_neutral_modelCul", ".pdf", sep = "")
ggsave(FileName, p1,width = 12,height = 4)
FileName <- paste(phypath,"1_neutral_modelCul", ".png", sep = "")
ggsave(FileName, p1,width = 12,height = 4)



#Phylogenetic signal #--
source("./phylo_Micro\\phyloSignal_and_phySigplot.R")
envphy = env
row.names(env) = env$ID
env$ID = NULL
# # No output from this function is saved to the specified path, because this calculation is a waste of time #---------
phypath2 = paste(phypath,"/phyloSignal/",sep = "")
dir.create(phypath)
phyloSignal(ps = psphy,
            Group = "group",
            Env = env,
            Path = phypath2)

Result = phySigPlot(ps = ps,group = "group ",env = env,path = phypath2)
#
# extract image
p2 = result[[1]] + mytheme1
p2
#- Extract plotting data
data = result[[2]]
head(data)

FileName<- paste(phypath,"2_phySigPlot", ".pdf", sep = "")
ggsave(FileName, p2,width = 15,height = 6)
FileNamepaste(phypath,"2_phySigPlot", ".csv", sep = "")
write.csv(data,FileName)



#-- Calculate the zero model #-------
source("./phylo_Micro\\nullModel1.R")

result <- nullModel(ps = psphy,
                    group="Group",
                    dist.method =  "bray",
                    gamma.method = "total",
                    transfer = "none",
                    null.model = "ecosphere"
)

#-- Group zero model run results
nullModeltabresult[[1]]

# proportion
ratiotabresult[[2]]
#- Statistics statistics difference
aovtabresult[[3]]

FileNamepaste(phypath,"3_nullModeltab", ".csv", sep = "")
write.csv(nullModeltab,FileName)

FileNamepaste(phypath,"3_ratiotab", ".csv", sep = "")
write.csv(ratiotab,FileName)

# FileNamepaste(phypath,"3_aovtab", ".csv", sep = "")
# write.csv(aovtab,FileName)




#--BNTI#----
source(".phylo_Micro\\bNTICul.R")

Result = bNTICul(ps = psphy,group = "group ",num = 999,thread = 1)

bNTI = result[[1]]
head(bNTI)


filename = paste(phypath,"/4_bNTI.csv",sep = "")
write.csv(bNTI, filename)


# - calculating RCbray# -- -- -- -- -- -- -- -- -- -- -
source("./phylo_Micro\\RCbary.R")

Result = RCbary(ps = psphy,num = 99,thread = 1)

RCbary = result[[1]]
head(RCbary)
filename = paste(phypath,"/5_RCb.csv",sep = "")
write.csv(RCbary,filename)

#--BetaNTI and RCbray collaborated on this image #---------
phypath = "./result_and_plot/16S_env_phylo_processing/Phylogenetic_analyse_spacies/"
source("./phylo_Micro\\bNTIRCPlot.R")

bNTI = read.csv(paste(phypath,"/4_bNTI.csv",sep = ""),row.names = 1)
head(bNTI)
# RCbray data read in, modify column name
RCb = read.csv(paste(phypath,"/5_RCb.csv",sep = ""),row.names = 1) %&gt; %
  Dplyr ::mutate(Sample_1 = Site2, Sample_2 = Site1)
head(RCb)

Result = bNTIRCPlot(ps = psphy,RCb = RCb,bNTI = bNTI,group = "group ")

# - bNTI out of the picture
p3result[[1]] + mytheme1
p3

# RCbary visualization
p4result[[2]] + mytheme1
p4
# Combine pictures BNTI, RCbray
p5 <- result[[3]]
p5
plotdata = result[[4]]
head(plotdata)

filename = paste(phypath,"/6_bNTI_RCbray.csv",sep = "")
write.csv(plotdata,filename)

FileName <- paste(phypath,"6_bNTI", ".pdf", sep = "")
ggsave(FileName, p3,width =8,height = 6)

FileName <- paste(phypath,"6_RCbary", ".pdf", sep = "")
ggsave(FileName, p4,width = 6,height = 6)

FileName <- paste(phypath,"6_BNTI_RCbray", ".pdf", sep = "")
ggsave(FileName, p5,width = 12,height = 8)

FileName <- paste(phypath,"6_bNTI", ".png", sep = "")
ggsave(FileName, p3,width =8,height = 6)

FileName <- paste(phypath,"6_RCbary", ".png", sep = "")
ggsave(FileName, p4,width = 6,height = 6)

FileName <- paste(phypath,"6_BNTI_RCbray", ".png", sep = "")
ggsave(FileName, p5,width = 12,height = 8)

#-- Environmental factors related to BetaNTI #---------
source("./phylo_Micro\\EnvCorbNTI.R")

#- Import bNTI functions
bNTIRC = read.csv(paste(phypath,"/6_bNTI_RCbray.csv",sep = ""),row.names = 1)
head(bNTIRC)

map = sample_data(psphy)
head(map)
plot = EnvCorbNTI(ps = psphy,
                  bNTIRC = bNTIRC,
                  group  = "Group",
                  env = envRDA
)

