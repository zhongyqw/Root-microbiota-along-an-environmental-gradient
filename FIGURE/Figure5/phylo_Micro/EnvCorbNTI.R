#
#
# library(tidyverse)
# library(vegan)
# library(picante)
# library(minpack.lm)
# library(FSA)
# library(eulerr)
# library(ggplot2)
# library(grid)
# require(Hmisc)#
# require(stats4)
#
#
# ps = readRDS("./ps.rds")
# ps
#
#
# #-导入bNTI函数
# bNTIRC = read.csv("../pipeline//bNTI_RCbray.csv",row.names = 1)
# head(bNTIRC)
#
# env = read.delim("../ori_data/env.txt",row.names = 1)
# head(env)
#
# plot = EnvCorbNTI(ps = ps,bNTIRC = bNTIRC,group  = "Group")
#
# ## 提取相关分析结果，总图
# plot[[1]]
# #提取单个
# # plot[[2]][1]

EnvCorbNTI = function(otu = NULL,
                      tax = NULL,
                      map = NULL,
                      tree = NULL,
                      ps = NULL,
                      bNTIRC = RCbNTI,
                      env = env,
                      group  = "Group"){

  ps = inputMicro(otu,tax,map,tree,ps,group  = group)
  ps
  #------------定义相关性分析函数
  # df = data
  Sams.mantel.test = function(df, seed=NULL) {
    # Run mantel test to see if there is a correlation
    delta.mat = df %>%
      select(Sample_1, Sample_2, delta) %>%
      spread(Sample_2, delta)
    rownames(delta.mat) = delta.mat$Sample_1
    delta.mat$Sample_1 = NULL
    delta.mat = delta.mat[names(sort(rowSums(!is.na(delta.mat)), decreasing = F)), names(sort(colSums(!is.na(delta.mat)), decreasing = T))]
    delta.mat = as.dist(delta.mat)

    bNTI.mat = df %>%
      select(Sample_1, Sample_2, bNTI) %>%
      spread(Sample_2, bNTI)
    rownames(bNTI.mat) = bNTI.mat$Sample_1
    bNTI.mat$Sample_1 = NULL
    bNTI.mat = bNTI.mat[names(sort(rowSums(!is.na(bNTI.mat)), decreasing = F)), names(sort(colSums(!is.na(bNTI.mat)), decreasing = T))]
    bNTI.mat = as.dist(bNTI.mat)
    if (!(is.null(seed))){
      set.seed(seed)
    }
    bNTI.mat[is.na(bNTI.mat)] = 0
    mantel.res = vegan::mantel(delta.mat, bNTI.mat)
    return(mantel.res)
  }




  set.seed(72)  # setting seed for reproducibility
  psrare = rarefy_even_depth(ps)
  #检查序列数量
  sample_sums(psrare)
  # 标准化数据
  ps.norm = transform_sample_counts(psrare, function(x) x/sum(x))




  map = as.data.frame(sample_data(psrare))
  # map = data.frame(row.names = map$id,id = map$id,Group = map$Group)
  mapE =merge(map,env,by = "row.names",all= FALSE)
  row.names(mapE) = mapE$Row.names
  mapE$Row.names = NULL
  mapE$ID = row.names(mapE)
  head(mapE)

  #---------合并环境变量数据
  # i = "Altitude..m."
  plot = list()
  for (i in colnames(env)) {

    colnames(mapE) = gsub(i,"XX",colnames(mapE))

    # Add in pH metadata
    pH.meta1=mapE %>%
      dplyr::select(ID, XX) %>%
      dplyr::rename(Sample_1 = ID, env1_1 = XX)

    pH.meta2= mapE%>%
      dplyr::select(ID, XX) %>%
      dplyr::rename(Sample_2 = ID, env1_2 = XX)

    data = dplyr::inner_join(bNTIRC, pH.meta1) %>%
      dplyr::inner_join(pH.meta2) %>%
      dplyr::mutate(delta = abs(env1_1-env1_2),
             crosstype = ifelse(Group_1 == Group_2, as.character(Group_1), "across"))
    head(data)
    data$crosstype
    # Run mantel test to see if there is a correlation
    pH.mantel = Sams.mantel.test(data, seed=72)
    head(data)
    # Plot
    p = ggplot(data, aes(x=delta, y=bNTI)) +
      geom_point(pch = 21) +
      # scale_shape_manual(values=LandUse.shapes) +
      geom_hline(yintercept = 2, linetype=2) +
      geom_hline(yintercept = -2, linetype=2) +
      # annotate("text", x=3.25, y=12.5, label=paste("r= ", round(pH.mantel$statistic, 3), "\n", "p= ", round(pH.mantel$signif, 3), sep="")) +
      labs(x=paste("",i), y="βNTI",title = paste("r= ", round(pH.mantel$statistic, 3), "p= ", round(pH.mantel$signif, 3))) +
      theme(legend.position = "none") +theme_bw()

    p
    plot[[i]] = p
    colnames(mapE) = gsub("XX",i,colnames(mapE))
  }

  library(ggpubr)
  p  = ggarrange(plotlist = plot, common.legend = TRUE, legend="right")
  p

  return(list(p,plot))
}

Sams.mantel.test = function(df, seed=NULL) {
  # Run mantel test to see if there is a correlation
  delta.mat = df %>%
    select(Sample_1, Sample_2, delta) %>%
    spread(Sample_2, delta)
  rownames(delta.mat) = delta.mat$Sample_1
  delta.mat$Sample_1 = NULL
  delta.mat = delta.mat[names(sort(rowSums(!is.na(delta.mat)), decreasing = F)), names(sort(colSums(!is.na(delta.mat)), decreasing = T))]
  delta.mat = as.dist(delta.mat)
  
  bNTI.mat = df %>%
    select(Sample_1, Sample_2, bNTI) %>%
    spread(Sample_2, bNTI)
  rownames(bNTI.mat) = bNTI.mat$Sample_1
  bNTI.mat$Sample_1 = NULL
  bNTI.mat = bNTI.mat[names(sort(rowSums(!is.na(bNTI.mat)), decreasing = F)), names(sort(colSums(!is.na(bNTI.mat)), decreasing = T))]
  bNTI.mat = as.dist(bNTI.mat)
  if (!(is.null(seed))){
    set.seed(seed)
  }
  bNTI.mat[is.na(bNTI.mat)] = 0
  mantel.res = vegan::mantel(delta.mat, bNTI.mat)
  return(mantel.res)
}
