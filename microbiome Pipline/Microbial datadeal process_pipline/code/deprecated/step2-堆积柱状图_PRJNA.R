library(pacman)
pacman::p_load(
  reshape2, ggplot2, RColorBrewer, yyplot, openxlsx,
  stringr, paletteer, trend, Hmisc
)

# 设置工作目录和参数
work_dir = "D:/Project/gutDBase/bracken_summary_GSE"
host = "human"
taxonomy = "phylum"
project = "GSE123649"
group_condition = "sex"

# 读取丰度和元数据
abundance = read.xlsx(
  file.path(
    work_dir, host, paste0(project, "_summary.xlsx")
  ),
  sheet = taxonomy,
  rowNames = TRUE
)

meta_data = read.xlsx(
  file.path(
    work_dir, host, paste0(project, "_summary.xlsx")
  ),
  sheet = "metaData",
  rowNames = TRUE
)

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
melt_data = melt_data[order(melt_data$Bacteroidota), ] # 按某一菌排序
melt_data$SampleID = rownames(melt_data)

# 转为长格式用于绘图
plot_data = melt(
  melt_data,
  id.vars = c("SampleID", colnames(meta_data)),
  measure.vars = phylumn_names,
  variable.name = "Phylum",
  value.name = "Abundance"
)

plot_data$SampleID = factor(plot_data$SampleID,
  levels = unique(plot_data$SampleID)
)
plot_data = plot_data[order(plot_data$Abundance), ]
plot_data$Abundance = as.numeric(plot_data$Abundance)

# 绘图准备
n_colors = length(unique(plot_data$Phylum))
mycolors = colorRampPalette(brewer.pal(8, "Set3"))(n_colors)
facet_var = as.formula(paste("~", group_condition))

# 绘图
gg3 = ggplot(plot_data, aes(x = SampleID)) +
  geom_bar(aes(y = Abundance, fill = Phylum), stat = "identity", width = 0.85) +
  facet_grid(facet_var, scales = "free_x", space = "free") +
  labs(x = "", y = "Relative Abundance") +
  guides(fill = guide_legend(reverse = FALSE)) +
  theme(
    text = element_text(family = "serif"),
    panel.background = element_blank(),
    panel.spacing = unit(0.1, "lines"),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    strip.text.x = element_text(size = 10, face = "bold", family = "serif"),
    strip.background.x = element_rect(fill = "#F6F6F6", colour = "#E0E0E0"),
    axis.title.y = element_text(size = 14, face = "bold"),
    legend.key.size = unit(15, "pt"),
    legend.title = element_text(size = 12, face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  scale_fill_manual(name = capitalize(taxonomy), values = mycolors)

gg3 = set_font(gg3, family = "serif")
