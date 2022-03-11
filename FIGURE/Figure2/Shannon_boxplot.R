library(ggpubr)
library(reshape2)
library(ggsci)

data=read.delim("Shannon.txt",header = T, sep="\t")

# Add the p-value between different groups
my_comparisons <- list(c("Endosp.", "Rhizopl."), c("Endosp.", "Rhizosp."), c("Endosp.", "Soil"),
                       c("Rhizopl.","Rhizosp."),c("Rhizopl.","Soil"),c("Rhizosp.","Soil"))
#Comparison between different groups
#label here means select significance label (asterisk)
p <- ggboxplot(data, x="part", y="shannon", color = "part",add = "jitter")+
  stat_compare_means(comparisons = my_comparisons,label = "p.signif")+  
  stat_compare_means(label.y = 3.5)+scale_y_continuous(limits = c(2.5,15))+
  theme_bw(base_size = 10)+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(axis.text.x = element_text(angle = 30, hjust = 1))+
  scale_color_aaas()
p
ggsave("part-shannon-boxplot.pdf",width=12,height=9,units="cm")

