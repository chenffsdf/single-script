#############################################################
# 脚本名称: 03.annotation.R
# 功能描述: 生成注释表格，原始表格+lookup注释表格（gene_type列补充）+GO/KEGG注释
#          （补充新的转录本GO注释）+events事件信息+统计不同种类event的数量，绘制
#           柱形图+加TPM+把gene symbol放到前面
# 作者: guowenbin & chenjiawei
# 创建日期: 2025-05-16
# 最后修改: 2025-05-16 (修改人: chenjiawei)
# 版本: 1.0
# 输入文件: *testing statistics.csv或Significant * statistics.csv;rnaseq_transfeat.csv;
#           interproscan_ann_nopw.tsv;AS_final.csv;intermediate_data.RData;lookup.RData;
#           GO_KEGG_term.csv
# 输出文件: *testing statistics.csv或Significant * statistics.csv;* events count.csv;
#           *events count.pdf(png)
# 依赖包: stringr;fs;dplyr;tidyr;GO.db;ggplot2
#############################################################

setwd("D:/华智种谷智创/4.项目/01.20250324_老鼠/20250515_小鼠项目后续分析")
library("stringr")
library("dplyr")
library("fs")
library("tidyr")
library("GO.db")
library("ggplot2")

sig_gene <- dir_ls(
  path = "./",
  regexp = "20250515_小鼠项目后续分析/.*/.*gene.* testing statistics.csv",
  recurse = TRUE,
  type = "file"
)
sig_transcript <- dir_ls(
  path = "./",
  regexp = "20250515_小鼠项目后续分析/.*/.*transcripts testing statistics.csv",
  recurse = TRUE,
  type = "file"
)


for (f in sig_gene){
  ###lookup注释表格（gene_type列补充）
  load('../lookup.RData')
  lookup <- as.data.frame(lookup)
  lookup[lookup == "NA"] <- NA
  map <- lookup[,c("gene_id", "transcript_id")]
  map <- unique(map)
  
  transfeat <- read.csv("../rnaseq_transfeat.csv")
  transfeat <- transfeat[,c("Transcript_ID","Coding_potentiality")]
  transfeat$Coding_potentiality[transfeat$Coding_potentiality == "Coding"] <- "protein_coding"
  names(transfeat)[1] = "transcript_id"
  lookup <- lookup %>%
    left_join(transfeat, by = "transcript_id") %>%
    mutate(gene_type = ifelse(is.na(gene_type), Coding_potentiality, gene_type)) %>%
    dplyr::select(-Coding_potentiality)
  dat_split1 <- split(lookup,lookup$gene_id)
  
  dat_split1 <- lapply(dat_split1,function(x){
    data.frame(t(apply(x, 2, function(y) paste(na.omit(unique(y[y!=''])),collapse = ';'))))
  })
  lookup <- do.call(rbind,dat_split1)
  df <- read.csv(f)
  #需更改
  anno <- df %>%
    left_join(lookup, by = c("target" = "gene_id"))
  term <- read.csv("GO_KEGG_term.csv")
  
  anno <- anno %>%
    left_join(term,by=c("target"="ID"))
  ###GO/KEGG注释（补充新的转录本GO注释
  interpro <- read.csv("../interproscan_ann_nopw.tsv",header = FALSE, sep = "\t")
  interpro <- interpro[,c("V1","V9","V14")]
  interpro <- interpro[interpro$V14 != "-",]
  interpro <- interpro %>%
    left_join(map,by=c("V1" = "transcript_id"))
  
  interpro <- interpro %>%
    group_by(gene_id) %>%               # 按 V1 列分组
    slice_min(V9, n = 1) %>%       # 保留每组中 V9 最小的行
    ungroup()                      # 取消分组
  interpro$V14 <- gsub("\\(InterPro\\)", "", interpro$V14)
  interpro$V14 <- gsub("\\(PANTHER\\)", "", interpro$V14)
  interpro <- tidyr::separate_rows(data = interpro, "V14",sep = "\\|")
  # 获取GO号对应的通路描述
  go_terms <- select(
    GO.db,
    keys = unique(interpro$V14),
    columns = c("GOID", "TERM"),  # 提取ID、描述和类别（BP/MF/CC）
    keytype = "GOID"
  )
  interpro <- interpro %>%
    left_join(go_terms, by = c("V14" = "GOID"))
  interpro$GO1 <- paste(interpro$V14, interpro$TERM, sep = "~")
  
  dat_split <- split(interpro,interpro$gene_id)
  dat_split <- lapply(dat_split,function(x){
    data.frame(t(apply(x, 2, function(y) paste(na.omit(unique(y[y!=''])),collapse = ';'))))
  })
  interpro <- do.call(rbind,dat_split)
  interpro <- interpro[,c("gene_id","GO1")]
  anno <- anno %>%
    left_join(interpro, by = c("target" = "gene_id")) %>%
    mutate(GO_term = ifelse(is.na(GO_term) | anno$GO_term == "", GO1, GO_term)) %>%
    dplyr::select(-GO1)
  
  
  ###添加events事件信息###
  ###统计不同种类event的数量，绘制柱形图###
  event <- read.csv("../AS_final.csv")
  #需要修改
  event <- event[,c("gene_id","event")]
  event <- unique(event)
  event_p <- event
  event_p$event <- sub(":.*", "", event_p$event)
  event_p <- anno %>%
    left_join(event_p, by=c("target"="gene_id"))
  
  event_counts <- data.frame(
    case = names(table(event_p$event)),  # 提取事件类型
    count = as.numeric(table(event_p$event))  # 提取频数
  )
  new_f <- sub("list and statistics\\.csv$", "events count.csv", f)
  #write.csv(event_counts, new_f, row.names = FALSE)
  pdf_name <- sub("/result/", "/figure/", 
                  sub("list and statistics\\.csv$", "events count.pdf", f))
  png_name <- sub("/result/", "/figure/", 
                  sub("list and statistics\\.csv$", "events count.png", f))
  
  
  p <- ggplot(event_counts, aes(x = case, y = count)) +
    geom_col(fill = "steelblue", width = 0.7) +
    geom_text(aes(label = count), vjust = -0.5, size = 3.5) +    # 添加数值标签
    labs(
      x = "Event Type", 
      y = "Count", 
      title = "Events count"
    ) +
    theme_classic() +
    theme(
      # 白色背景和黑色边框
      panel.background = element_rect(fill = "white", color = "black"),
      plot.background = element_rect(fill = "white", color = "black"),
      # 坐标轴刻度标签大小
      axis.text.x = element_text(size = 12, colour = "black"),  # x轴刻度大小
      axis.text.y = element_text(size = 12, colour = "black"),             # y轴刻度大小
      panel.grid.minor = element_blank(),
      # 坐标轴标题大小
      axis.title.x = element_text(size = 16, colour = "black"),  # x轴标题大小
      axis.title.y = element_text(size = 16, colour = "black"),  # y轴标题大小
      plot.title = element_text(hjust = 0.5, face = "bold",size = 20)
    )
  p
  #ggsave(pdf_name,p,width = 7,height = 5)
  #ggsave(png_name,p,width = 7,height = 5)
  # 需要修改
  dat_split <- split(event,event$gene_id)
  
  dat_split <- lapply(dat_split,function(x){
    data.frame(t(apply(x, 2, function(y) paste(na.omit(unique(y[y!=''])),collapse = ';'))))
  })
  
  event <- do.call(rbind,dat_split)
  # 需要修改
  anno <- anno %>%
    left_join(event, by=c("target"="gene_id"))
  
  #### 添加TPM
  part <- strsplit(f, "/")[[1]][2]
  if (part == "FBMSC_MBMSC"){
    intermediate <- "../20250506_3D_结果1/20250507_mouse_3D_结果1/data/intermediate_data.RData"
  } else if (part == "FGH_FBMSC") {
    intermediate <- "../20250506_3D_结果2/20250507_mouse_3D_结果2/data/intermediate_data.RData"
  }
  load(intermediate)
  genes_TPM <- intermediate_data$genes_TPM
  # 提取target id
  genes <- anno$target
  # 只取跟id匹配的行
  tpm <- genes_TPM[genes,]
  colnames(tpm) <- paste0('TPM:',colnames(tpm))
  tpm <- as.data.frame(tpm)
  tpm$id <- rownames(tpm)
  anno <- anno %>%
    left_join(tpm, by=c("target"="id"))
  anno <- anno %>% 
    dplyr::select(1, transcript_id, gene_name, everything())
  anno$event <- str_trunc(
    anno$event,
    width = 32750,  # 最大长度
    side = "right",  # 从右侧截断
    ellipsis = "······"  # 截断后添加的符号
  )
  write.csv(anno,f,row.names = FALSE,na = '')
}

for (f in sig_transcript){
  load('../lookup.RData')
  lookup <- as.data.frame(lookup)
  lookup[lookup == "NA"] <- NA
  
  ###GO/KEGG注释（补充新的转录本GO注释
  #需更改
  df <- read.csv(f)
  anno <- df %>%
    left_join(lookup, by = c("target" = "transcript_id"))
  term <- read.csv("GO_KEGG_term.csv")
  anno <- anno %>%
    left_join(term, by=c("gene_id"="ID"))

  interpro <- read.csv("../interproscan_ann_nopw.tsv",header = FALSE, sep = "\t")
  interpro <- interpro[,c("V1","V9","V14")]
  interpro <- interpro[interpro$V14 != "-",]
  interpro <- interpro %>%
    left_join(map,by=c("V1" = "transcript_id"))
  
  interpro <- interpro %>%
    group_by(gene_id) %>%               # 按 V1 列分组
    slice_min(V9, n = 1) %>%       # 保留每组中 V9 最小的行
    ungroup()                      # 取消分组
  interpro$V14 <- gsub("\\(InterPro\\)", "", interpro$V14)
  interpro$V14 <- gsub("\\(PANTHER\\)", "", interpro$V14)
  interpro <- tidyr::separate_rows(data = interpro, "V14",sep = "\\|")
  # 获取GO号对应的通路描述
  go_terms <- select(
    GO.db,
    keys = unique(interpro$V14),
    columns = c("GOID", "TERM"),  # 提取ID、描述和类别（BP/MF/CC）
    keytype = "GOID"
  )
  interpro <- interpro %>%
    left_join(go_terms, by = c("V14" = "GOID"))
  interpro$GO1 <- paste(interpro$V14, interpro$TERM, sep = "~")
  
  dat_split <- split(interpro,interpro$gene_id)
  dat_split <- lapply(dat_split,function(x){
    data.frame(t(apply(x, 2, function(y) paste(na.omit(unique(y[y!=''])),collapse = ';'))))
  })
  interpro <- do.call(rbind,dat_split)
  interpro <- interpro[,c("gene_id","GO1")]
  anno <- anno %>%
    left_join(interpro, by = "gene_id") %>%
    mutate(GO_term = ifelse(is.na(GO_term) | anno$GO_term == "", GO1, GO_term)) %>%
    dplyr::select(-GO1)
  
  
  ###添加events事件信息###
  ###统计不同种类event的数量，绘制柱形图###
  event <- read.csv("../AS_final.csv")
  #需要修改
  event <- event[,c("transcript_id","event")]
  event <- unique(event)
  event_p <- event
  event_p$event <- sub(":.*", "", event_p$event)
  event_p <- anno %>%
    left_join(event_p, by=c("target"="transcript_id"))
  
  event_counts <- data.frame(
    case = names(table(event_p$event)),  # 提取事件类型
    count = as.numeric(table(event_p$event))  # 提取频数
  )
  new_f <- sub("list and statistics\\.csv$", "events count.csv", f)
  #write.csv(event_counts, new_f, row.names = FALSE)
  pdf_name <- sub("/result/", "/figure/", 
                  sub("list and statistics\\.csv$", "events count.pdf", f))
  png_name <- sub("/result/", "/figure/", 
                  sub("list and statistics\\.csv$", "events count.png", f))
  
  
  p <- ggplot(event_counts, aes(x = case, y = count)) +
    geom_col(fill = "steelblue", width = 0.7) +
    geom_text(aes(label = count), vjust = -0.5, size = 3.5) +    # 添加数值标签
    labs(
      x = "Event Type", 
      y = "Count", 
      title = "Events count"
    ) +
    theme_classic() +
    theme(
      # 白色背景和黑色边框
      panel.background = element_rect(fill = "white", color = "black"),
      plot.background = element_rect(fill = "white", color = "black"),
      # 坐标轴刻度标签大小
      axis.text.x = element_text(size = 12, colour = "black"),  # x轴刻度大小
      axis.text.y = element_text(size = 12, colour = "black"),             # y轴刻度大小
      panel.grid.minor = element_blank(),
      # 坐标轴标题大小
      axis.title.x = element_text(size = 16, colour = "black"),  # x轴标题大小
      axis.title.y = element_text(size = 16, colour = "black"),  # y轴标题大小
      plot.title = element_text(hjust = 0.5, face = "bold",size = 20)
    )
  p
  #ggsave(pdf_name,p,width = 7,height = 5)
  #ggsave(png_name,p,width = 7,height = 5)
  # 需要修改
  dat_split <- split(event,event$transcript_id)
  
  dat_split <- lapply(dat_split,function(x){
    data.frame(t(apply(x, 2, function(y) paste(na.omit(unique(y[y!=''])),collapse = ';'))))
  })
  
  event <- do.call(rbind,dat_split)
  # 需要修改
  anno <- anno %>%
    left_join(event, by=c("target"="transcript_id"))
  
  ###添加TPM矩阵
  part <- strsplit(f, "/")[[1]][2]
  if (part == "FBMSC_MBMSC"){
    intermediate <- "../20250506_3D_结果1/20250507_mouse_3D_结果1/data/intermediate_data.RData"
  } else if (part == "FGH_FBMSC") {
    intermediate <- "../20250506_3D_结果2/20250507_mouse_3D_结果2/data/intermediate_data.RData"
  }
  load(intermediate)
  trans_TPM <- intermediate_data$trans_TPM
  # 提取target id
  trans <- anno$target
  # 只取跟id匹配的行
  tpm <- trans_TPM[trans,]
  colnames(tpm) <- paste0('TPM:',colnames(tpm))
  tpm <- as.data.frame(tpm)
  tpm$id <- rownames(tpm)
  anno <- anno %>%
    left_join(tpm, by=c("target"="id"))
  anno <- anno %>% 
    dplyr::select(1, gene_id, gene_name, everything())
  anno$event <- str_trunc(
    anno$event,
    width = 32750,  # 最大长度
    side = "right",  # 从右侧截断
    ellipsis = "······"  # 截断后添加的符号
  )
  write.csv(anno,f,row.names = FALSE,na = '')
  
}














