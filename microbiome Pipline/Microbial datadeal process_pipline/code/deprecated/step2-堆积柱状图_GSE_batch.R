pacman::p_load(
  reshape2, ggplot2, RColorBrewer, yyplot, openxlsx,
  stringr, paletteer, trend, Hmisc
)

plot_microbiome_barplot = function(
    accession,
    group_info,
    taxonomy = "phylum",
    work_dir,
    meta_file,
    host,
    order_by = "Bacteroidota",
    font_family = "serif") {
  pie_meta = read.csv(meta_file, stringsAsFactors = FALSE)
  sub = subset(pie_meta, accession == accession & info == group_info)
  group_list = split(sub$sample, sub$value)

  if (!file.exists(file.path(
    work_dir,
    host, paste0(accession, "_summary.xlsx")
  ))) {
    cat(paste(
      "跳过：", accession, "-", group_info,
      "没有对应的 summary 文件"
    ))
    return(NULL)
  }

  # 读取丰度和元数据
  abundance = read.xlsx(
    file.path(work_dir, host, paste0(accession, "_summary.xlsx")),
    sheet = taxonomy,
    rowNames = TRUE
  )

  meta_data = read.xlsx(
    file.path(work_dir, host, paste0(accession, "_summary.xlsx")),
    sheet = "metaData",
    rowNames = TRUE
  )

  if (length(meta_data$Sample.Name) != ncol(abundance)) {
    cat(paste(
      "跳过：", accession, "-", group_info,
      "下载数据不全"
    ))
    return(NULL)
  }

  # 数据归一化
  abundance_norm = scale(abundance, center = FALSE, scale = colSums(abundance))
  rownames(abundance_norm) = str_split(
    rownames(abundance_norm),
    paste0(strsplit(taxonomy, "")[[1]][1], "__"),
    simplify = TRUE
  )[, 2]

  phylumn_names = rownames(abundance_norm)

  # 合并数据
  melt_data = cbind(meta_data, t(abundance_norm))
  if (!(order_by %in% rownames(abundance_norm))) {
    cat(paste(
      "跳过：", accession, "-", group_info,
      "中未检测到", order_by, "，无法排序"
    ))
    return(NULL)
  }
  melt_data = melt_data[order(melt_data[[order_by]]), ]
  melt_data$SampleID = rownames(melt_data)

  # 长格式
  plot_data = melt(
    melt_data,
    id.vars = c("SampleID", colnames(meta_data)),
    measure.vars = phylumn_names,
    variable.name = "Phylum",
    value.name = "Abundance"
  )

  plot_data$Sample.Name = factor(plot_data$Sample.Name,
    levels = unique(plot_data$Sample.Name)
  )
  plot_data = plot_data[order(plot_data$Abundance), ]
  plot_data$Abundance = as.numeric(plot_data$Abundance)

  # 添加 GSM 分组
  plot_data$GSMgroup = NA_character_
  for (grp in names(group_list)) {
    plot_data$GSMgroup[plot_data$Sample.Name %in% group_list[[grp]]] = grp
  }
  plot_data = subset(plot_data, !is.na(GSMgroup))

  # 设置颜色
  n_colors = length(unique(plot_data$Phylum))
  mycolors = colorRampPalette(brewer.pal(8, "Set3"))(n_colors)

  facet_var = as.formula("~ GSMgroup")

  gg3 = ggplot(plot_data, aes(x = Sample.Name)) +
    geom_bar(aes(y = Abundance, fill = Phylum), stat = "identity", width = 0.85) +
    facet_grid(facet_var, scales = "free_x", space = "free") +
    labs(x = "", y = "Relative Abundance") +
    guides(fill = guide_legend(reverse = FALSE)) +
    theme(
      text = element_text(family = font_family),
      panel.background = element_blank(),
      panel.spacing = unit(0.1, "lines"),
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank(),
      strip.text.x = element_text(size = 10, face = "bold", family = font_family),
      strip.background.x = element_rect(fill = "#F6F6F6", colour = "#E0E0E0"),
      axis.title.y = element_text(size = 14, face = "bold"),
      legend.key.size = unit(15, "pt"),
      legend.title = element_text(size = 12, face = "bold"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    scale_fill_manual(name = capitalize(taxonomy), values = mycolors)

  gg3 = set_font(gg3, family = font_family)
  return(gg3)
}

# 设置工作目录和参数
work_dir = "D:/Project/gutDBase/bracken_summary_GSE"
host = "human"
taxonomy = "phylum"

# 读取完整 meta 文件
pie_meta = read.csv("D:/Project/gutDBase/metadata/hum_pie.csv",
  stringsAsFactors = FALSE
)

# 获取唯一的 accession 和 info 组合
combo_table = unique(pie_meta[, c("accession", "info")])

# 创建输出文件夹（如果不存在）
output_dir = file.path("D:/Project/gutDBase/Abundance_Barplot", host)
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# 遍历组合
for (i in seq_len(nrow(combo_table))) {
  acc = combo_table$accession[i]
  info = combo_table$info[i]
  cat("绘图中：", acc, "-", info, "\n")

  try({
    plot_result = plot_microbiome_barplot(
      accession = acc,
      group_info = info,
      taxonomy = "phylum",
      work_dir = work_dir,
      meta_file = "D:/Project/gutDBase/metadata/hum_pie.csv",
      host = host,
      order_by = "Bacteroidota",
      font_family = "serif"
    )
    file_name = file.path(output_dir, paste0(acc, "_", info, ".png"))
    if (!is.null(plot_result)) {
      ggsave(file_name, plot_result, width = 12, height = 6, dpi = 300)
    }
  })
}
