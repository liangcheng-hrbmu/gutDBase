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

# 设置路径与读取数据
work_dir = "D:/Project/gutDBase/bracken_summary_GSE"
host = "human"
project = "GSE123649"
taxonomy = "species"
group_condition = "sex"

abundance = read.xlsx(file.path(work_dir, host, paste0(project, "_summary.xlsx")), sheet = taxonomy, rowNames = TRUE)
meta_data = read.xlsx(file.path(work_dir, host, paste0(project, "_summary.xlsx")), sheet = "metaData", rowNames = TRUE)

# 数据标准化
abundance_norm = scale(abundance, center = FALSE, scale = colSums(abundance))
rownames(abundance_norm) = str_split(rownames(abundance_norm),
  paste0(strsplit(taxonomy, "")[[1]][1], "__"),
  simplify = TRUE
)[, 2]

sample_id = rownames(meta_data)
group_by_info = meta_data[[group_condition]]

# 多样性指数计算
chao_ace = estimateR(t(abundance))
chao1 = chao_ace[2, ]
shannon = diversity(t(abundance_norm), index = "shannon")
simpson = diversity(t(abundance_norm), index = "simpson")

# 整理绘图数据
data_plot = data.frame(sample_id, group_by_info, chao1, shannon, simpson)
colnames(data_plot)[2] = group_condition
data_plot[[group_condition]] = factor(data_plot[[group_condition]], levels = unique(data_plot[[group_condition]]))

# α-多样性可视化函数
plot_alpha_diversity = function(plot_data, diversity_type, groupby, mycolor) {
  compairs = combn(levels(plot_data[, groupby]), 2)
  comparisons = list()

  for (i in 1:ncol(compairs)) {
    wt_a = plot_data[, diversity_type][plot_data[, groupby] == compairs[1, i]]
    wt_b = plot_data[, diversity_type][plot_data[, groupby] == compairs[2, i]]
    test_res = wilcox.test(wt_a, wt_b, paired = FALSE, exact = FALSE)
    if (test_res$p.value < 0.05) {
      comparisons = c(comparisons, list(compairs[, i]))
    }
  }

  ggplot(plot_data, aes(plot_data[, groupby], plot_data[, diversity_type])) +
    geom_flat_violin(
      aes(color = plot_data[, groupby], fill = plot_data[, groupby]),
      position = position_nudge(x = .25), alpha = 0.3, linewidth = 0.3
    ) +
    geom_jitter(aes(color = plot_data[, groupby]), width = 0.1, size = 1.1) +
    geom_boxplot(
      aes(color = plot_data[, groupby]),
      width = .1, position = position_nudge(x = 0.25), size = 0.3
    ) +
    labs(title = diversity_type, x = groupby, y = diversity_type) +
    theme_bw() +
    theme(
      text = element_text(family = "serif"),
      aspect.ratio = 1,
      panel.border = element_rect(fill = NA, linewidth = 0.3),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.title = element_text(face = "bold", size = 12),
      legend.title = element_text(face = "bold")
    ) +
    geom_signif(
      comparisons = comparisons, map_signif_level = TRUE,
      test = wilcox.test, textsize = 4, size = 0.2,
      vjust = 0.6, step_increase = 0.05, color = "#3D3D3D"
    ) +
    scale_fill_manual(name = str_to_title(gsub("_", " ", groupby)), values = mycolor) +
    scale_color_manual(name = str_to_title(gsub("_", " ", groupby)), values = mycolor)
}

# 设置颜色
mycolor = pal_npg("nrc")(10)[1:7]

# 多图合并
alpha_agent =
  plot_alpha_diversity(data_plot, "chao1", group_condition, mycolor) +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(hjust = 0.5, vjust = 0.5), plot.title = element_blank(), legend.position = "none") +
  plot_alpha_diversity(data_plot, "shannon", group_condition, mycolor) +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(hjust = 0.5, vjust = 0.5), plot.title = element_blank(), legend.position = "none") +
  plot_alpha_diversity(data_plot, "simpson", group_condition, mycolor) +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(hjust = 0.5, vjust = 0.5), plot.title = element_blank()) +
  plot_annotation(title = "α-diversity", theme = theme(plot.title = element_text(hjust = 0.5)))

alpha_agent
