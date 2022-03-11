library(ggplot2) 
library(ggbiplot)
library(ggsci)

pcdata <-  read.delim("../data/pcoa.txt",header = T, sep="\t")

#Add confidence circle
p = ggplot(pcdata,aes(PC1,PC2,color=part))+geom_point(size=2.5)+ 
    theme_bw(base_size=10)+xlab("PCo1 16.57%")+ylab("PCo2 5.49%")+
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
    geom_hline(yintercept=0,linetype="dashed",size=0.3)+
    geom_vline(xintercept=0,linetype="dashed",size=0.3)+
    stat_ellipse(level = 0.95)+scale_color_aaas()
p
ggsave("PCoA_part_BC.pdf",width=12,height=8,units="cm")


