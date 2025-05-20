#############################################################
# 脚本名称: Structure_Abundance_PS_plot.R
# 功能描述: 生成指定基因的结构图、Abundance图和PS图，并整合成1张图
# 作者: guowenbin & chenjiawei
# 创建日期: 2025-05-20
# 最后修改: 2025-05-20 (修改人: chenjiawei)
# 版本: 1.0
# 输入文件: eTrans;intermediate_data.RData;txi_trans.RData;mouse_rtd_merged_final_transfix
# 输出文件: Abundance gene.png(pdf);PS gene.png(pdf);gene_structure.png(pdf);gene_all.png(pdf)
# 依赖包: eTrans;cowplot;rtracklayer;ggplot2;ThreeDRNAseq
#############################################################

setwd("D:/华智种谷智创/4.项目/01.20250324_老鼠/20250519_小鼠项目后续分析")
sourceDir <- function(path, trace = TRUE, ...) {
  for (nm in list.files(path, pattern = '*.R')) {
    #if(trace) cat(nm,":")
    source(file.path(path, nm), ...)
    #if(trace) cat("/n")
  }
}
sourceDir('eTrans/R')
library("cowplot")
library("rtracklayer")
library("ggplot2")
library("ThreeDRNAseq")


projs <- c("FBMSC_MBMSC","FGH_FBMSC")

for (proj in projs){
  proj <- "FGH_FBMSC"
  if (proj == "FBMSC_MBMSC"){
    load("../20250506_3D_结果1/20250507_mouse_3D_结果1/data/intermediate_data.RData")
    load("../20250506_3D_结果1/20250507_mouse_3D_结果1/data/txi_trans.RData")
  } else if (proj == "FGH_FBMSC"){
    load("../20250506_3D_结果2/20250507_mouse_3D_结果2/data/intermediate_data.RData")
    load("../20250506_3D_结果2/20250507_mouse_3D_结果2/data/txi_trans.RData") 
  }
  DTU_trans <- intermediate_data$DTU_trans
  DTU_trans <- DTU_trans[DTU_trans$adj.pval<0.05,]
  # DTU_trans <- DTU_trans[order(DTU_trans$adj.pval,decreasing = F),]
  DTU_trans <- DTU_trans[grep(pattern = 'chr',DTU_trans$target),]
  genes <- unique(gsub('[.].*','',DTU_trans$target))
  print(genes)
  
  gtf <- import('../mouse_rtd_merged_final_transfix.gtf')
  for (gene in genes){
    gene <- 'chr12G012664'
    gene_folder <- paste0("20250519_小鼠项目后续分析/",proj,"/",gene)
    if (!dir.exists(gene_folder)) {
      dir.create(gene_folder, recursive = TRUE, showWarnings = FALSE)
    }
    gr <- gtf[gtf$gene_id==gene]
    rtracklayer::export(gr,paste0(gene_folder,"/",'gene_',gene,'.gtf'))
    trans_num <- length(unique(as.data.frame(gr)$transcript_id))
    p1 <- plotTrans(gr = gr,exon_color = "black",intron_color = 'black',exon_width = 5,intron_width = 1,
                    showCDS = T,cds_width_ratio = 1.8,
                    x_ticks = 10,arrow_ticks = 30)
    p1
    
    ggsave(paste0(gene_folder,"/",gene,"_structure.png"), 
           plot = p1, width = 24, height = 20)
    ggsave(paste0(gene_folder,"/",gene,"_structure.pdf"), 
           plot = p1, width = 24, height = 20)
    
    
    
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
    
    # 保留全部的转录本
    g.pr <- plotAbundance(data.exp = txi_trans$abundance,
                          gene = gene,
                          mapping = mapping,
                          genes.ann = genes.ann,
                          trans.ann = trans.ann,
                          trans.expressed = NULL,
                          reps = metatable$label,
                          y.lab = 'TPM')
    g.pr
    ggsave(paste0(gene_folder,"/Abundance ",gene,".png"),
           plot = g.pr, width = 18, height = 12, dpi = 300)
    ggsave(paste0(gene_folder,"/Abundance ",gene,".pdf"), 
           plot = g.pr, width = 18, height = 12)
    
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
    ggsave(paste0(gene_folder,"/PS ",gene,".png"),
           plot = g.ps, width = 18, height = 12, dpi = 300)
    ggsave(paste0(gene_folder,"/PS ",gene,".pdf"), 
           plot = g.ps, width = 18, height = 12)
    
    
    p <- plot_grid(
      p1,
      plot_grid(g.pr, g.ps, ncol = 2),
      nrow = 2,
      rel_heights = c(1, 1)
    )
    p
    ggsave(paste0(gene_folder,"/",gene,"_all.png"),
           plot = p, width = 32, height = 24, dpi = 300)
    ggsave(paste0(gene_folder,"/",gene,"_all.pdf"), 
           plot = p, width = 32, height = 24)
  }
}

