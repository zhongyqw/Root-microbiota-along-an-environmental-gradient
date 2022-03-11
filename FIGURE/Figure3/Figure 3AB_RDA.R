
install.packages("devtools")
library(devtools)
devtools::install_github("gavinsimpson/ggvegan")
library(ggvegan)
library(vegan)

#read data
otu.tab <- read.table("../data/RB.core.asv.16S.txt",header=T,row.names=1,sep="\t")
env.data <- read.table("../data/RB.env2.txt",header=T,fill=T,row.names=1,sep="\t")

#transform data
otu <- t(otu.tab)
#data normolization
env.data <- env.data[3:16]
env.data.log <- log1p(env.data)##
env <- na.omit(env.data.log)

###hellinger transform
otu.hell <- decostand(otu, "hellinger")
# Extract only those ID in common between the two tables
idx = rownames(otu.hell) %in% rownames(env)
env = env[idx,]
env = env[rownames(otu.hell),]
env = na.omit(env)
#DCA analysis
sel <- decorana(otu.hell)
sel

db_rda <- capscale(otu.hell~., env, distance = 'bray', add = TRUE)
db_rda

db_rda.scaling1 <- summary(db_rda, scaling = 1)
db_rda.scaling1

#RsquareAdj() Extract  R2，
r2 <- RsquareAdj(db_rda)
db_rda_noadj <- r2$r.squared #The original R2
db_rda_noadj
db_rda_adj <- r2$adj.r.squared       #After correction R2
db_rda_adj

#Restraint shaft correction
db_rda_exp_adj <- db_rda_adj * db_rda$CCA$eig/sum(db_rda$CCA$eig)
db_rda_exp_adj #The first two axes are explanatory rates
db_rda_eig_adj <- db_rda_exp_adj * db_rda$tot.chi
db_rda_eig_adj

#Displacement test
#The permutation test for all constrained axes, i.e., the global test, is based on 999 permutations
db_rda_test <- anova.cca(db_rda, permutations = 999)
db_rda_test
#Each constraint axis was tested one by one based on 999 permutations
db_rda_test_axis <- anova.cca(db_rda, by = 'axis', permutations = 999)
db_rda_test_axis
#Each environmental factor was tested one by one based on 999 substitutions
db_rda_test_term <- anova.cca(db_rda, by = 'term', permutations = 999)
db_rda_test_term
#P-value correction
db_rda_test_axis$`Pr(>F)` <- p.adjust(db_rda_test_axis$`Pr(>F)`, method = 'bonferroni')
db_rda_test_axis$`Pr(>F)`

#############
#forward choose
db_rda_forward_pr <- ordiR2step(capscale(otu.hell~1, env, distance = 'bray', add = TRUE), scope = formula(db_rda), R2scope = TRUE, direction = 'forward', permutations = 999)
#
db_rda_forward_pr <- capscale(otu.hell~MAP+NO3.+pH+SOC+TP, env, distance = 'bray', add = TRUE)
db_rda_forward_pr

#View details
summary(db_rda_forward_pr, scaling = 1)
#Compare the difference of R2 after correction before and after selection
RsquareAdj(db_rda)$adj.r.squared
RsquareAdj(db_rda_forward_pr)$adj.r.squared

#Global check of all constraint axes, 999 permutations
db_rda_forward_pr_test <- anova.cca(db_rda_forward_pr, permutations = 999)
db_rda_forward_pr_test
#Each constraint axis was tested one by one, and 999 permutations were performed
db_rda_forward_pr_test_axis <- anova.cca(db_rda_forward_pr, by = 'axis', permutations = 999)
db_rda_forward_pr_test_axis
#Each constraint axis was tested one by one, and 999 permutations were performed
db_rda_forward_pr_test_term <- anova.cca(db_rda_forward_pr, by = 'term', permutations = 999)
db_rda_forward_pr_test_term

##########

#Quadrat scores were calculated using species weights and
db_rda_forward_pr_site.scaling1 <- scores(db_rda_forward_pr, choices = 1:2, scaling = 1, display = 'wa')
#write.table(data.frame(db_rda_forward_pr_site.scaling1), 'db_rda_forward_pr_site.scaling1.txt', sep = '\t', col.names = NA, quote = FALSE)
#Species variable (response variable) score
db_rda_forward_pr_sp.scaling1 <- scores(db_rda_forward_pr, choices = 1:2, scaling = 1, display = 'sp')
#write.table(data.frame(db_rda_forward_pr_sp.scaling1), 'db_rda_forward_pr_sp.scaling1.txt', sep = '\t', col.names = NA, quote = FALSE)
#Environmental variables (explanatory variables) score
db_rda_forward_pr_env.scaling1 <- scores(db_rda_forward_pr, choices = 1:2, scaling = 1, display = 'bp')
#write.table(data.frame(db_rda_forward_pr_env.scaling1), 'db_rda_forward_pr_env.scaling1.txt', sep = '\t', col.names = NA, quote = FALSE)

######

db_rda_vp <- varpart(otu.hell,  ~ MAT+MAP, ~SOC+TN+TP+NO3.+NH4.+pH+BD,~AGB+BGB+R.OC
                     +R.N, data=env, scale = FALSE, add = TRUE)
db_rda_vp

plot(db_rda_vp, digits = 3, Xnames = c('X1', 'X2','X3'), bg = c('blue', 'red','green'))

##################################
##Extract RDA data
smry <- summary(db_rda_forward_pr, scaling = 1)
df1  <- data.frame(smry$sites[,1:2])       # RDA1 and RDA2
df2  <- data.frame(smry$species[,1:2])     # loadings for RDA1 and RDA2
df3  <- data.frame(smry$biplot[,1:2])
#Plus grouping information
group<-read.table("group.txt",header = T, row.names = 1,sep = "\t")
sub_group = group[group[,1]=="Rhizopl.",] 
df11<-merge(df1,sub_group,by = "row.names", all = TRUE)

#The constraint-axis interpretation rate was calculated after R2 correction
exp_adj <- RsquareAdj(db_rda_forward_pr)$adj.r.squared * db_rda_forward_pr$CCA$eig/sum(db_rda_forward_pr$CCA$eig)
rda1_exp <- paste('RDA1:', round(exp_adj[1]*100, 2), '%')
rda2_exp <- paste('RDA1:', round(exp_adj[2]*100, 2), '%')

library(ggrepel) 
library(ggplot2)
#Painting reproduction
rda.plot <- ggplot(df11, aes(x=CAP1, y=CAP2,colour=type)) + geom_point(size=2.5,alpha=0.8)+
  geom_hline(yintercept=0, linetype="dotted") +
  geom_vline(xintercept=0, linetype="dotted")
#Add the arrow
rda.biplot <- rda.plot +
  geom_segment(data=df3, aes(x=0, xend=CAP1, y=0, yend=CAP2),color="blue", arrow=arrow(length=unit(0.03,"npc")))+
  geom_text_repel(data=df3, aes(x=CAP1,y=CAP2,label=rownames(df3)),color="blue", size=3)+
  labs(x = rda1_exp, y = rda2_exp) +
  theme_bw(base_size=10)+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
rda.biplot

library(ggsci)
rda.biplot+scale_color_aaas()
ggsave("soil.core.rda.pdf",width=15,height=9.5,units="cm")

