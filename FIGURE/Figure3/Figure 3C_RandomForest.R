library(randomForest)
data=read.delim("RB.env2.txt",header = T, sep="\t")
set.seed(315)
rf = randomForest(data=data, PC1~MAT+MAP+
                    pH+SOC+TN+TP+NH4.+NO3.+R.OC+R.N+AGB+BGB+BD, 
                  importance=TRUE, proximity=TRUE, ntree = 1000,na.action=na.omit) 
rf = randomForest(data=data, PC2~MAT+MAP+
                    pH+SOC+TN+TP+NH4.+NO3.+R.OC+R.N+AGB+BGB+BD,importance=TRUE, proximity=TRUE, ntree = 1000,na.action=na.omit) 
print(rf)
#The importance of the independent variable
importance(rf)
#Draw an order of importance for variables
varImpPlot(rf)
#feature importance
imp= as.data.frame(rf$importance)
imp = imp[order(imp[,1],decreasing = T),] 
head(imp)
write.table(imp,file = "imp_S_PC1.txt",quote = F,sep = '\t', row.names = T, col.names = T)
write.table(imp,file = "imp_S_PC2.txt",quote = F,sep = '\t', row.names = T, col.names = T)

library(ggplot2)
names(imp) <- c("IncMSE", "Inc")
imp$X <- row.names(imp)
imp$X=factor(imp$X,levels = imp$X)
pc1=ggplot(data = imp, mapping = aes(x=X,y=IncMSE))+
  geom_bar(stat="identity",fill="#E7B800")+
  theme_bw(base_size=8)+xlab("factors")+ylab("Importance (Increase in MSE %)")+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
pc1
ggsave("S.PC1.pdf",width=10,height=6,units="cm")