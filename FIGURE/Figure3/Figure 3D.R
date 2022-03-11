library(ggplot2)
library(ggpubr)
library(reshape2)
library(ggsci)
library(ggpmisc)

data=read.table("clipboard",header = T, sep="\t")
#data<-read.table(pipe("pbpaste"),header = T, sep="\t")

data_long=melt(data, id.vars = c("MAP","type"),
variable.names= "part",
value.name="dist")
head(data_long)              
data_long=na.omit(data_long)

###fit
formula <- y ~ x
p<-ggplot(data_long, aes(x = MAP, y = dist))+geom_point(alpha=0.6)+
  geom_smooth(method = "lm", aes(group=1),formula = formula) +
  facet_wrap(variable~.,nrow=2,scales = "free")+
  #stat_cor(method = "pearson") +
  stat_poly_eq(
    aes(label =  paste(..p.value.label.., ..adj.rr.label.., sep = "~~~~")),
    formula = formula, parse = TRUE,eq.with.lhs = F,
  )+
  theme_bw(base_size = 8)+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
p
ggsave("distance with MAP.pdf",width=16,height=8,units="cm")
#
formula <- y ~ x
p<-ggplot(data_long, aes(x = MAP, y = dist,colour=type))+geom_point(alpha=0.6)+
  geom_smooth(method = "lm", aes(group=1),formula = formula) +
  facet_wrap(variable~.,nrow=2,scales = "free")+
  #stat_cor(method = "pearson") +
  theme_bw(base_size = 8)+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  scale_color_aaas()
p
ggsave("distance with MAP.pdf",width=16,height=8,units="cm")
