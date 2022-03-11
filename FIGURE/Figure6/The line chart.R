library(ggplot2)
library(ggpubr)
library(reshape2)
library(ggsci)


data=read.delim("../data/Network stability.txt",header = T, sep="\t")
data_long=melt(data, id.vars = c("part","index"),
               variable.names= "part",
               value.name="value")
head(data_long)              
data_long=na.omit(data_long)
data_long$index<-factor(data_long$index,levels =c("num.vertices","modularity","ratio of negative:positive cohesion"))
p<-ggplot(data_long, aes(x = variable, y = value,group=1,colour=variable))+
  geom_line(linetype = "dashed",color="red")+geom_point()+
  facet_grid(index~part,scales = "free")+
  theme_bw(base_size = 8)+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  scale_color_npg()
p
ggsave("network.pdf",width=20,height=10,units="cm")
