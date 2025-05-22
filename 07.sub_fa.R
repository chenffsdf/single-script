#############################################################
# 脚本名称: 07.sub_fa.R
# 功能描述: 生成指定基因的转录本fa序列
# 作者: chenjiawei
# 创建日期: 2025-05-22
# 最后修改: 2025-05-22 (修改人: chenjiawei)http://127.0.0.1:35369/graphics/plot_zoom_png?width=1170&height=824
# 版本: 1.0
# 输入文件: mouse_rtd_merged_final.fasta;intermediate_data.RData
# 输出文件: gene.fa
# 依赖包: Biostrings
#############################################################

setwd("D:/华智种谷智创/4.项目/01.20250324_老鼠/20250520_小鼠项目_提取转录本序列")
library("Biostrings")

load("../20250506_3D_结果2/20250507_mouse_3D_结果2/data/intermediate_data.RData")
map <- intermediate_data$mapping

genes <- c("chr2G001630","chr2G002076","chr4G004356")


for (gene in genes){
  # gene <- "chr2G001630"
  trans <- map[map$GENEID == gene,"TXNAME"]
  fa <- "../mouse_rtd_merged_final.fasta"
  RNA <- readDNAStringSet(fa)
  rna_sub <- RNA[trans]
  Biostrings::writeXStringSet(rna_sub,paste0("20250520_小鼠项目_提取转录本序列/",gene,"/",gene,".fa"))
}

