
# dir.amp(ps0)
# theme.amp(theme = T)
# package.amp()
# color.amp()
#-扩增子分析需要导入的R包
package.amp <- function(){
  library(phyloseq)
  library(tidyverse)
  library(ggClusterNet)
  library(EasyStat)
  library(fs)
  library(ggthemes)
  library(RColorBrewer)#调色板调用包
  library(magrittr)
  library(MicrobiotaProcess)
  library(ggsignif)
  library(ggtree)
  library(ggtreeExtra)
  library(ggstar)
  library(MicrobiotaProcess)
  library(ggnewscale)
  library(grid)
}

dir.amp <- function(ps0){
  # # 建立结构保存一级目录#--------
  result_path <- paste("./","/result_and_plot/",sep = "")
  fs::dir_create(result_path)
  #---构建结果保存文件夹#---------
  tax.1 = c("Fungi",
    "fungi",
    "K__Fungi",
    "k__Fungi",
    "d__Fungi",
    "d__fungi",
    "d__Eukaryota",
    "Eukaryota",
    "K:Fungi",
    "k:Fungi",
    "d:Fungi",
    "d:fungi",
    "d:Eukaryota",
    "Eukaryota",
    "D:Eukaryota"
    )
  TFdir.f <- as.data.frame(table(
    phyloseq::tax_table(ps0)[,1]))[,2][as.data.frame(table(phyloseq::tax_table(ps0)[,1]))[,1] %in%
                               tax.1] > 10
  if (length(TFdir.f) != 0) {
    print("ITS")
    res1path <- paste(result_path,"/Base_diversity_ITS",sep = "")
    id = tax.1
  }
  
  tax.2 = c("Bacteria",
            "K__Bacteria",
            "k__Bacteria",
            "d__Bacteria",
            "k:Bacteria",
            "bacteria",
            "D:bacteria",
            "D__bacteria",
            "d__bacteria",
            "d:prokaryotes",
            "d__prokaryotes",
            "prokaryotes"
            )
  
  TFdir.b <- as.data.frame(table(phyloseq::tax_table(ps0)[,1]))[,2][as.data.frame(table(phyloseq::tax_table(ps0)[,1]))[,1] %in%
                                                            tax.2 ] > 10
  
  if (length(TFdir.b) != 0) {
    print("16s")
    res1path <- paste(result_path,"/Base_diversity_16s",sep = "")
    id = tax.2
  }
  
  fs::dir_create(res1path)
  return(list(res1path,id))
}

#-包括几个常用的主题



theme_my = function(ps = ps0) {
  mytheme1 = ggplot2::theme_bw() + ggplot2::theme(
    panel.background=  ggplot2::element_blank(),
    panel.grid=  ggplot2::element_blank(),
    legend.position="right",
    legend.title =  ggplot2::element_blank(),
    legend.background= ggplot2::element_blank(),
    legend.key= ggplot2::element_blank(),
    plot.title =  ggplot2::element_text(vjust = -8.5,hjust = 0.1),
    axis.title.y =  ggplot2::element_text(size = 24,face = "bold",colour = "black"),
    axis.title.x = ggplot2::element_text(size = 24,face = "bold",colour = "black"),
    axis.text =  ggplot2::element_text(size = 20,face = "bold"),
    axis.text.x =  ggplot2::element_text(colour = "black",size = 14),
    axis.text.y =  ggplot2::element_text(colour = "black",size = 14),
    legend.text =  ggplot2::element_text(size = 15,face = "bold")
  )
  
  mytheme2 = ggplot2::theme_bw() + ggplot2::theme(
    panel.background=  ggplot2::element_blank(),
    panel.grid=  ggplot2::element_blank(),
    legend.position="right",
    
    legend.title =  ggplot2::element_blank(),
    legend.background=  ggplot2::element_blank(),
    legend.key= ggplot2::element_blank(),
    plot.title =  ggplot2::element_text(vjust = -8.5,hjust = 0.1),
    axis.title.y =  ggplot2::element_text(size = 24,face = "bold",colour = "black"),
    axis.title.x = ggplot2::element_text(size = 24,face = "bold",colour = "black"),
    axis.text =  ggplot2::element_text(size = 20,face = "bold"),
    axis.text.x =  ggplot2::element_text(colour = "black",size = 14,angle = 90),
    axis.text.y =  ggplot2::element_text(colour = "black",size = 14),
    legend.text =  ggplot2::element_text(size = 15,face = "bold")
  )
  
  #--设定颜色
  gnum <- unique(phyloseq::sample_data(ps)$Group) %>% length()
  
  # scales::show_col(RColorBrewer::brewer.pal(9,"Set1"))
    if (gnum < 10 ) {
      colset1 <- RColorBrewer::brewer.pal(9,"Set1")
    } else {
      colset1 <- colorRampPalette(RColorBrewer::brewer.pal(9,"Set1"))(gnum)
    }
    #设定颜色#------------

    colset2 <- RColorBrewer::brewer.pal(12,"Paired")
    colset3 <- c(RColorBrewer::brewer.pal(11,"Set1"),RColorBrewer::brewer.pal(9,"Pastel1"))
    colset4 = colset3
  

  return(list(mytheme1,mytheme2,colset1,colset2,colset3,colset4))
  
}





