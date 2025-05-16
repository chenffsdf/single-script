#############################################################
# 脚本名称: 01.make_gene_gtf.R
# 功能描述: 生成指定基因的gtf文件
# 作者: guowenbin & chenjiawei
# 创建日期: 2025-05-16
# 最后修改: 2025-05-16 (修改人: chenjiawei)
# 版本: 1.0
# 输入文件: mouse_rtd_merged_final_transfix.gtf
# 输出文件: gene.gtf
# 依赖包: rtracklayer
#############################################################

library(rtracklayer)

gtf <- import('mouse_rtd_merged_final_transfix.gtf')
genes1 <- c("Igkv5-39",
            "Igkv10-96",
            "Igkv19-93",
            "Ighg2b")
for (gene in genes1){
  gr <- gtf[gtf$gene_name %in% gene]
  gr <- gr[gr$type %in% c('CDS','exon')]
  folder_path <- paste0("20250515_小鼠项目后续分析/FBMSC_MBMSC/", gene)
  if (!dir.exists(folder_path)) {
    dir.create(folder_path, recursive = TRUE, showWarnings = FALSE)
  }
  rtracklayer::export(gr,paste0("20250515_小鼠项目后续分析/FBMSC_MBMSC/",gene,"/",'gene_',gene,'.gtf'))
}