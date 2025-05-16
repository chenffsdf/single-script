#############################################################
# 脚本名称: 02.Abundance_PS_plot.R
# 功能描述: 生成指定基因的Abundance图和PS图
# 作者: guowenbin & chenjiawei
# 创建日期: 2025-05-16
# 最后修改: 2025-05-16 (修改人: chenjiawei)
# 版本: 1.0
# 输入文件: intermediate_data.RData;txi_trans.RData
# 输出文件: Abundance gene.png(pdf);PS gene.png(pdf)
# 依赖包: ThreeDRNAseq;ggplot2
#############################################################

setwd("D:/华智种谷智创/4.项目/01.20250324_老鼠/20250515_小鼠项目后续分析")
library("ggplot2")
library("ThreeDRNAseq")

genes1 <- c("Igkv5-39",
            "Igkv10-96",
            "Igkv19-93",
            "Ighg2b")
genes1_es <- c("ENSMUSG00000076569",
               "ENSMUSG00000094420",
               "ENSMUSG00000098814",
               "ENSMUSG00000076613")


for (i in 1:length(genes1)){
  load("../20250506_3D_结果1/20250507_mouse_3D_结果1/data/intermediate_data.RData")
  load("../20250506_3D_结果1/20250507_mouse_3D_结果1/data/txi_trans.RData")
  
  DE_genes <- intermediate_data$DE_genes
  DAS_genes <- intermediate_data$DAS_genes
  DE_trans <- intermediate_data$DE_trans
  DTU_trans <- intermediate_data$DTU_trans
  target_high <- intermediate_data$target_high
  metatable <- intermediate_data$samples
  metatable$label <- metatable$condition
  mapping <- intermediate_data$mapping
  
  DE.genes <- unique(DE_genes$target)
  DAS.genes <- unique(DAS_genes$target)
  DE.trans <- unique(DE_trans$target)
  DTU.trans <- unique(DTU_trans$target)
  
  genes.ann <- set2(DE.genes,DAS.genes)
  names(genes.ann) <- c('DEonly','DE&DAS','DASonly')
  genes.ann <- plyr::ldply(genes.ann,cbind)[,c(2,1)]
  colnames(genes.ann) <- c('target','annotation')
  
  
  trans.ann <- set2(DE.trans,DTU.trans)
  names(trans.ann) <- c('DEonly','DE&DTU','DTUonly')
  trans.ann <- plyr::ldply(trans.ann,cbind)[,c(2,1)]
  colnames(trans.ann) <- c('target','annotation')
  
  ################################################################################
  ##give a gene name
  gene <- genes1_es[i]
  gene_name <- genes1[i]
  
  ##----->> Abundance plot
  # #只保留高表达
  # g.pr <- plotAbundance(data.exp = txi_trans$abundance[target_high$trans_high,],
  #                       gene = gene,
  #                       mapping = mapping[target_high$trans_high,],
  #                       genes.ann = genes.ann,
  #                       trans.ann = trans.ann,
  #                       trans.expressed = NULL,
  #                       reps = metatable$label,
  #                       y.lab = 'TPM')+
  #   labs(title = paste0(gene,' (',gene_name,')'))
  # g.pr
  
  # 保留全部的转录本
  g.pr <- plotAbundance(data.exp = txi_trans$abundance,
                        gene = gene,
                        mapping = mapping,
                        genes.ann = genes.ann,
                        trans.ann = trans.ann,
                        trans.expressed = NULL,
                        reps = metatable$label,
                        y.lab = 'TPM')+
    labs(title = paste0(gene,' (',gene_name,')'))
  g.pr
  ggsave(paste0("20250515_小鼠项目后续分析/FBMSC_MBMSC/",gene_name,"/","Abundance ",gene,".png"),
         plot = g.pr, width = 8, height = 6, dpi = 300)
  ggsave(paste0("20250515_小鼠项目后续分析/FBMSC_MBMSC/",gene_name,"/","Abundance ",gene,".pdf"), 
         plot = g.pr, width = 8, height = 6)
  
  ##----->> PS plot
  g.ps <- plotPS(data.exp = txi_trans$abundance[target_high$trans_high,],
                 gene = gene,
                 mapping = mapping[target_high$trans_high,],
                 genes.ann = genes.ann,
                 trans.ann = trans.ann,
                 trans.expressed = NULL,
                 reps = metatable$label,
                 y.lab = 'PS')
  g.ps
  ggsave(paste0("20250515_小鼠项目后续分析/FBMSC_MBMSC/",gene_name,"/","PS ",gene,".png"),
         plot = g.ps, width = 8, height = 6, dpi = 300)
  ggsave(paste0("20250515_小鼠项目后续分析/FBMSC_MBMSC/",gene_name,"/","PS ",gene,".pdf"), 
         plot = g.ps, width = 8, height = 6)
}
