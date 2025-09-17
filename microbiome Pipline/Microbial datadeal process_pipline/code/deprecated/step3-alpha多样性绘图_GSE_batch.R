# 加载必要的 R 包
library(pacman)
pacman::p_load(
  vegan, ggplot2, ggpubr, rstatix, ggsci, ggsignif,
  reshape2, gridExtra, scales, cowplot, dplyr, readr, patchwork,
  stringr, lavaan, smooth, Hmisc, openxlsx, readr
)

# 定义小提琴图函数
"%||%" = function(a, b) {
  if (!is.null(a)) a else b
}

geom_flat_violin = function(mapping = NULL, data = NULL, stat = "ydensity",
                            position = "dodge", trim = TRUE, scale = "area",
                            show.legend = NA, inherit.aes = TRUE, ...) {
  layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomFlatViolin,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(trim = trim, scale = scale, ...)
  )
}

GeomFlatViolin = ggproto(
  "GeomFlatViolin", Geom,
  setup_data = function(data, params) {
    data$width = data$width %||% params$width %||% (resolution(data$x, FALSE) * 0.9)
    data %>%
      group_by(group) %>%
      mutate(
        ymin = min(y),
        ymax = max(y),
        xmin = x,
        xmax = x + width / 2
      )
  },
  draw_group = function(data, panel_scales, coord) {
    data = transform(data,
      xminv = x,
      xmaxv = x + violinwidth * (xmax - x)
    )
    newdata = rbind(
      plyr::arrange(transform(data, x = xminv), y),
      plyr::arrange(transform(data, x = xmaxv), -y)
    )
    newdata = rbind(newdata, newdata[1, ])
    ggplot2:::ggname(
      "geom_flat_violin",
      GeomPolygon$draw_panel(newdata, panel_scales, coord)
    )
  },
  draw_key = draw_key_polygon,
  default_aes = aes(
    weight = 1, colour = "grey20", fill = "white",
    linewidth = 0.5, alpha = NA, linetype = "solid"
  ),
  required_aes = c("x", "y")
)

plot_alpha_diversity_all = function(
    accession,
    group_info,
    taxonomy,
    work_dir,
    host,
    meta_file,
    font_family = "serif") {
  pie_meta = read.csv(meta_file, stringsAsFactors = FALSE)
  sub = subset(pie_meta, accession == accession & info == group_info)
  if (nrow(sub) == 0) {
    cat("跳过：", accession, "-", group_info, "无匹配分组信息\n")
    return(NULL)
  }

  file_path = file.path(work_dir, host, paste0(accession, "_summary.xlsx"))
  if (!file.exists(file_path)) {
    cat("跳过：", accession, "-", group_info, "缺少 summary 文件\n")
    return(NULL)
  }

  abundance = read.xlsx(file_path, sheet = taxonomy, rowNames = TRUE)
  meta_data = read.xlsx(file_path, sheet = "metaData", rowNames = TRUE)

  # 判断样本数量是否一致
  if (length(meta_data$Sample.Name) != ncol(abundance)) {
    cat("跳过：", accession, "-", group_info, "下载数据不全\n")
    return(NULL)
  }

  # 数据标准化
  abundance_norm = scale(abundance, center = FALSE, scale = colSums(abundance))
  rownames(abundance_norm) = str_split(
    rownames(abundance_norm),
    paste0(strsplit(taxonomy, "")[[1]][1], "__"),
    simplify = TRUE
  )[, 2]

  # 提取所有组的 sample ID
  group_list = split(sub$sample, sub$value)
  sample_included = unlist(group_list)

  # 只保留指定分组样本
  abundance = abundance[, sample_included, drop = FALSE]
  abundance_norm = abundance_norm[, sample_included, drop = FALSE]
  meta_data = meta_data[sample_included, , drop = FALSE]

  # 构造分组因子
  group_factor = rep(NA_character_, length(sample_included))
  names(group_factor) = sample_included
  for (grp in names(group_list)) {
    group_factor[group_list[[grp]]] = grp
  }
  group_factor = factor(group_factor[sample_included])

  # 计算多样性
  chao_ace = estimateR(t(abundance))
  chao1 = chao_ace[2, sample_included]
  shannon = diversity(t(abundance_norm), index = "shannon")[sample_included]
  simpson = diversity(t(abundance_norm), index = "simpson")[sample_included]

  data_plot = data.frame(
    sample_id = sample_included,
    group = group_factor,
    chao1 = chao1,
    shannon = shannon,
    simpson = simpson
  )
  colnames(data_plot)[2] = group_info
  data_plot[[group_info]] = factor(data_plot[[group_info]])

  # 子图函数
  plot_alpha_diversity = function(plot_data, diversity_type, groupby) {
    compairs = combn(levels(plot_data[[groupby]]), 2)
    comparisons = list()
    for (i in 1:ncol(compairs)) {
      a = plot_data[[diversity_type]][plot_data[[groupby]] == compairs[1, i]]
      b = plot_data[[diversity_type]][plot_data[[groupby]] == compairs[2, i]]
      if (length(a) > 1 & length(b) > 1) {
        test_res = wilcox.test(a, b)
        if (test_res$p.value < 0.05) {
          comparisons = c(comparisons, list(compairs[, i]))
        }
      }
    }

    ggplot(plot_data, aes_string(groupby, diversity_type)) +
      geom_flat_violin(aes_string(fill = groupby), alpha = 0.3) +
      geom_jitter(aes_string(color = groupby), width = 0.1, size = 1.1) +
      geom_boxplot(aes_string(color = groupby), width = .1, size = 0.3) +
      labs(title = diversity_type, x = groupby, y = diversity_type) +
      theme_bw() +
      theme(
        text = element_text(family = font_family),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold")
      ) +
      geom_signif(
        comparisons = comparisons, map_signif_level = TRUE,
        test = wilcox.test, textsize = 4, size = 0.2,
        step_increase = 0.05
      ) +
      scale_fill_npg() +
      scale_color_npg()
  }

  # 组合输出
  alpha_plot =
    plot_alpha_diversity(data_plot, "chao1", group_info) +
    theme(legend.position = "none") +
    plot_alpha_diversity(data_plot, "shannon", group_info) +
    theme(legend.position = "none") +
    plot_alpha_diversity(data_plot, "simpson", group_info) +
    plot_annotation(title = "α-diversity", theme = theme(plot.title = element_text(hjust = 0.5)))

  return(alpha_plot)
}

# 设置路径
meta_file = "D:/Project/gutDBase/metadata/hum_pie.csv"
work_dir = "D:/Project/gutDBase/bracken_summary_GSE"
host = "human"
taxonomy = "species"
output_dir = file.path("D:/Project/gutDBase/Alpha_Diversity", host)
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# 读取元信息表，获取所有 accession-info 组合
pie_meta = read.csv(meta_file, stringsAsFactors = FALSE)
combo_table = unique(pie_meta[, c("accession", "info")])

# 遍历所有组合并绘图保存
for (i in seq_len(nrow(combo_table))) {
  acc = combo_table$accession[i]
  info = combo_table$info[i]
  cat("绘图中：", acc, "-", info, "\n")

  file_path = file.path(work_dir, host, paste0(acc, "_summary.xlsx"))
  if (!file.exists(file_path)) {
    cat("跳过：文件不存在 -", file_path, "\n")
    next
  }

  try({
    gg = plot_alpha_diversity_all(
      accession = acc,
      group_info = info,
      taxonomy = taxonomy,
      work_dir = work_dir,
      host = host,
      meta_file = meta_file
    )
    if (!is.null(gg)) {
      ggsave(
        filename = file.path(output_dir, paste0(acc, "_", info, "_alpha.png")),
        plot = gg, width = 12, height = 5, dpi = 300
      )
    }
  })
}
