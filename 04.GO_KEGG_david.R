#############################################################
# 脚本名称: 04.GO_KEGG_david.R
# 功能描述: 整理David(https://davidbioinformatics.nih.gov/)富集结果文件，绘制TOP20的富集条形图
# 作者: guowenbin & chenjiawei
# 创建日期: 2025-05-16
# 最后修改: 2025-05-16 (修改人: chenjiawei)
# 版本: 1.0
# 输入文件: *GO(KEGG) annotation simplified.csv
# 输出文件: *GO top 20.png(pdf)
# 依赖包: tidyverse
#############################################################

library("tidyverse")
setwd("D:/华智种谷智创/9.test/01.富集分析/01.网站david/20250328_3D_mouse-1")
p_or_q <- "p"

gene_cases <- c("DE genes","DAS genes")
trans_cases <- c("DE transcripts", "DTU transcripts")
cases <- c("DE genes","DAS genes","DE transcripts", "DTU transcripts")

##GO分析
for (case in cases){
  go <- read.csv(paste0("result/", case," GO annotation simplified.csv"))
  go <- go[go$Category != "Category", ] 
  go$Category <- sapply(strsplit(go$Category, "_"), function(x) {
    if (length(x) >= 2) x[2] else NA
  })
  go$PValue <- as.numeric(go$PValue)
  go$FDR <- as.numeric(go$FDR)
  go$neg_log10_FDR <- -log10(go$FDR)
  go$neg_log10_PValue <- -log10(go$PValue)
  write.csv(go,paste0("result/", case," GO annotation simplified.csv"),row.names = F)
  if (p_or_q == "p"){
    go <- go[go$PValue <= 0.05, ]
    if (nrow(go) > 20) {
      top20_go <- go[order(go$PValue), ][1:20, ]
    } else {
      top20_go <- go
    }
    if (nrow(top20_go) > 0) {
      top20_go$Term1 <- sapply(strsplit(top20_go$Term, "~"), function(x) {
        if (length(x) >= 2) x[2] else NA
      })
      p1 <- ggplot(top20_go, aes(x = neg_log10_PValue, y = reorder(Term1, neg_log10_PValue))) +
        geom_col(aes(fill = Category), width = 0.7) +  # 使用柱状图
        facet_grid(Category ~ ., scales = "free_y", space = "free_y") +  # 按Category分面
        labs(
          x = "-log10(PValue)",
          y = "GO Term",
          title = "TOP 20 GO enrichment",
        ) +
        theme_bw() +
        theme(
          # 设置全局字体为Arial
          text = element_text(family = "sans", color = "black"),
          # 主标题设置（加大字号）
          plot.title = element_text(size = 12, hjust = 0.5, color = "black"),  # hjust居中
          # 坐标轴标题设置（X轴和Y轴）
          axis.title.x = element_text(size = 12, color = "black"),  # 横坐标标题大小
          axis.title.y = element_text(size = 12, color = "black"),  # 纵坐标标题大小
          axis.text.x = element_text(color = "black"),
          axis.text.y = element_text(size = 10, color = "black"),
          strip.text.y = element_text(color = "black"),  # 分面标签水平显示
          panel.spacing = unit(0.5, "lines"),    # 调整分面间距
          legend.position = "none"
        ) +
        scale_fill_brewer(palette = "Set1")        # 使用彩色填充
      ggsave(paste0("figure/",case, " GO top 20.png"), p1, height = 8, width = 6)
      ggsave(paste0("figure/",case, " GO top 20.pdf"), p1, height = 8, width = 6)
    }
  } else if (p_or_q == "q"){
    go <- go[go$FDR <= 0.05, ]  
    if (nrow(go) > 20) {
      top20_go <- go[order(go$FDR), ][1:20, ]
    } else {
      top20_go <- go
    }
    if (nrow(top20_go) > 0){
      top20_go$Term1 <- sapply(strsplit(top20_go$Term, "~"), function(x) {
        if (length(x) >= 2) x[2] else NA
      })
      p1 <- ggplot(top20_go, aes(x = neg_log10_FDR, y = reorder(Term1, neg_log10_FDR))) +
        geom_col(aes(fill = Category), width = 0.7) +  # 使用柱状图
        facet_grid(Category ~ ., scales = "free_y", space = "free_y") +  # 按Category分面
        labs(
          x = "-log10(FDR)", 
          y = "GO Term",
          title = "TOP 20 GO enrichment",
        ) +
        theme_bw() +
        theme(
          # 设置全局字体为Arial
          text = element_text(family = "sans", color = "black"),
          # 主标题设置（加大字号）
          plot.title = element_text(size = 12, hjust = 0.5, color = "black"),  # hjust居中
          # 坐标轴标题设置（X轴和Y轴）
          axis.title.x = element_text(size = 12, color = "black"),  # 横坐标标题大小
          axis.title.y = element_text(size = 12, color = "black"),  # 纵坐标标题大小
          axis.text.x = element_text(color = "black"),
          axis.text.y = element_text(size = 10, color = "black"),
          strip.text.y = element_text(color = "black"),  # 分面标签水平显示
          panel.spacing = unit(0.5, "lines"),    # 调整分面间距
          legend.position = "none"
        ) +
        scale_fill_brewer(palette = "Set1")        # 使用彩色填充
      ggsave(paste0("figure/",case, " GO top 20.png"), p1, height = 8, width = 6)
      ggsave(paste0("figure/",case, " GO top 20.pdf"), p1, height = 8, width = 6)
    }
  } else {
    warning("p_or_q必须是'p'或'q'")
  }
}


##KEGG富集
for (case in cases){
  kegg <- read.csv(paste0("result/", case," KEGG annotation simplified.csv"))
  kegg$PValue <- as.numeric(kegg$PValue)
  kegg$FDR <- as.numeric(kegg$FDR)
  kegg$neg_log10_FDR <- -log10(kegg$FDR)
  kegg$neg_log10_PValue <- -log10(kegg$PValue)
  write.csv(kegg,paste0("result/", case," KEGG annotation simplified.csv"),row.names = F)
  if (p_or_q == "p"){
    kegg <- kegg[kegg$PValue <= 0.05, ]  
    if (nrow(kegg) > 20) {
      top20_kegg <- kegg[order(kegg$PValue), ][1:20, ]
    } else {
      top20_kegg <- kegg
    }
    if (nrow(top20_kegg) > 0) {
      top20_kegg$Term1 <- sapply(strsplit(top20_kegg$Term, ":"), function(x) {
        if (length(x) >= 2) x[2] else NA
      })
      p2 <- ggplot(top20_kegg, aes(x = neg_log10_PValue, 
                                   y = reorder(Term1, neg_log10_PValue))) +
        geom_col(width = 0.7, fill="red") +  # 使用柱状图
        labs(
          x = "-log10(PValue)",
          y = "KEGG Term",
          title = "TOP 20 KEGG enrichment",
        ) +
        theme_bw() +
        theme(
          # 设置全局字体为Arial
          text = element_text(family = "sans", color = "black"),
          # 主标题设置（加大字号）
          plot.title = element_text(size = 12, hjust = 0.5, color = "black"),  # hjust居中
          # 坐标轴标题设置（X轴和Y轴）
          axis.title.x = element_text(size = 12, color = "black"),  # 横坐标标题大小
          axis.title.y = element_text(size = 12, color = "black"),  # 纵坐标标题大小
          axis.text.x = element_text(color = "black"),
          axis.text.y = element_text(size = 10, color = "black"),
          legend.position = "none"
        )
      ggsave(paste0("figure/",case, " KEGG top 20.png"), p2, height = 8, width = 6)
      ggsave(paste0("figure/",case, " KEGG top 20.pdf"), p2, height = 8, width = 6)
    }
  } else if (p_or_q == "q"){
    kegg <- kegg[kegg$FDR <= 0.05, ]  
    if (nrow(kegg) > 20) {
      top20_kegg <- kegg[order(kegg$PValue), ][1:20, ]
    } else {
      top20_kegg <- kegg
    }
    if (nrow(top20_kegg) > 0) {
      top20_kegg$Term1 <- sapply(strsplit(top20_kegg$Term, ":"), function(x) {
        if (length(x) >= 2) x[2] else NA
      })
      p2 <- ggplot(top20_kegg, aes(x = neg_log10_FDR, 
                                   y = reorder(Term1, neg_log10_FDR))) +
        geom_col(width = 0.7, fill="red") +  # 使用柱状图
        labs(
          x = "-log10(FDR)", 
          y = "KEGG Term",
          title = "TOP 20 KEGG enrichment",
        ) +
        theme_bw() +
        theme(
          # 设置全局字体为Arial
          text = element_text(family = "sans", color = "black"),
          # 主标题设置（加大字号）
          plot.title = element_text(size = 12, hjust = 0.5, color = "black"),  # hjust居中
          # 坐标轴标题设置（X轴和Y轴）
          axis.title.x = element_text(size = 12, color = "black"),  # 横坐标标题大小
          axis.title.y = element_text(size = 12, color = "black"),  # 纵坐标标题大小
          axis.text.x = element_text(color = "black"),
          axis.text.y = element_text(size = 10, color = "black"),
          legend.position = "none"
        )
      ggsave(paste0("figure/",case, " KEGG top 20.png"), p2, height = 8, width = 6)
      ggsave(paste0("figure/",case, " KEGG top 20.pdf"), p2, height = 8, width = 6)
    }
  } else {
    warning("p_or_q必须是'p'或'q'")
  }
}
