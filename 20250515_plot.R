setwd("D:/华智种谷智创/4.项目/01.20250324_老鼠/20250515_小鼠项目后续分析")
library(ggplot2)
# 我按这个分析逻辑筛了下靶基因，麻烦先出下“FGH vs FBMSC”分析样下面6个基因的分析图：Fth1、Foxp1、Postn、fn1、Pparg、chr12G012664，谢谢。
# 过老师，“FBMSC vs MBMSC”麻烦出下这4个基因的图Igkv5-39、Igkv10-96、Igkv19-93、Ighg2b（只发生了DE事件，没发生DAS事件-做了DE事件与DAS事件的交集只有2个基因，这组分析数据没“FGH vs FBMSC”明显）

# Significant DE genes list and statistics”这个表，麻烦放宽下阈值：卡p值<0.05 & FC >1.5与FC<0.75 - 补充两个分析的4个总表，把总表发他就行

genes1 <- c("Igkv5-39",
            "Igkv10-96",
            "Igkv19-93",
            "Ighg2b")
genes1_es <- c("ENSMUSG00000076569",
               "ENSMUSG00000094420",
               "ENSMUSG00000098814",
               "ENSMUSG00000076613")
# genes2 <- c("Fth1","Foxp1","Postn","Fn1","Pparg","chr12G012664")
genes2 <- c("Fth1",
            "Foxp1",
            "Postn",
            "Fn1",
            "Pparg",
            "chr12G012664")
genes2_es <- c("ENSMUSG00000024661",
               "ENSMUSG00000030067",
               "ENSMUSG00000027750",
               "ENSMUSG00000026193",
               "ENSMUSG00000000440",
               "chr12G012664"
               )
## Plot example
library(rtracklayer)
gtf <- import('mouse_rtd_merged_final_transfix.gtf')

for (gene in genes1){
  gr <- gtf[gtf$gene_name %in% gene]
  gr <- gr[gr$type %in% c('CDS','exon')]
  
  ##高表达量的转录本
  # trans_high <- intermediate_data$target_high$trans_high
  # gr <- gr[gr$transcript_id %in% trans_high]
  folder_path <- paste0("20250515_小鼠项目后续分析/FBMSC_MBMSC/", gene)
  if (!dir.exists(folder_path)) {
    dir.create(folder_path, recursive = TRUE, showWarnings = FALSE)
  }
  rtracklayer::export(gr,paste0("20250515_小鼠项目后续分析/FBMSC_MBMSC/",gene,"/",'gene_',gene,'.gtf'))
}
for (gene in genes2){
  gr <- gtf[gtf$gene_name %in% gene]
  gr <- gr[gr$type %in% c('CDS','exon')]
  
  ##高表达量的转录本
  # trans_high <- intermediate_data$target_high$trans_high
  # gr <- gr[gr$transcript_id %in% trans_high]
  folder_path <- paste0("20250515_小鼠项目后续分析/FGH_FBMSC/", gene)
  if (!dir.exists(folder_path)) {
    dir.create(folder_path, recursive = TRUE, showWarnings = FALSE)
  }
  rtracklayer::export(gr,paste0("20250515_小鼠项目后续分析/FGH_FBMSC/",gene,"/",'gene_',gene,'.gtf'))
}



################################################################################
##----->> generate annotation of DE and/or DAS, DE and/or DTU
library("ThreeDRNAseq")

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
  ##----->> profile plot of TPM or read counts
  
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
  
  
  #都保留
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


for (i in 1:length(genes2)){
  i <- 6
  load("../20250506_3D_结果2/20250507_mouse_3D_结果2/data/intermediate_data.RData")
  load("../20250506_3D_结果2/20250507_mouse_3D_结果2/data/txi_trans.RData")
  
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
  gene <- genes2_es[i]
  gene_name <- genes2[i]
  ##----->> profile plot of TPM or read counts
  
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
  
  
  #都保留
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
  
  ggsave(paste0("20250515_小鼠项目后续分析/FGH_FBMSC/",gene_name,"/","Abundance ",gene,".png"),
         plot = g.pr, width = 16, height = 8, dpi = 300)
  ggsave(paste0("20250515_小鼠项目后续分析/FGH_FBMSC/",gene_name,"/","Abundance ",gene,".pdf"), 
         plot = g.pr, width = 16, height = 8)
  
  # png(paste0("Abundance ",gene,".png"), width = 2000, height = 1200, res = 300)
  # print(g.pr)
  # dev.off()
  # pdf(paste0("Abundance ",gene,".pdf"), width = 10, height = 6)
  # print(g.pr)
  # dev.off()
  
  ################################################################################
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
  ggsave(paste0("20250515_小鼠项目后续分析/FGH_FBMSC/",gene_name,"/","PS ",gene,".png"),
         plot = g.ps, width = 16, height = 8, dpi = 300)
  ggsave(paste0("20250515_小鼠项目后续分析/FGH_FBMSC/",gene_name,"/","PS ",gene,".pdf"), 
         plot = g.ps, width = 16, height = 8)
  
}

