# 加载相关软件包
library(pacman)
pacman::p_load(
  tidyverse,
  ggpubr,
  rstatix,
  ggsci,
  ggsignif,
  reshape2,
  ggplot2,
  vegan,
  ggpubr,
  ape,
  patchwork,
  RColorBrewer,
  openxlsx
)

##### 数据读取
# 读取微生物属级别丰度数据
setwd(
  "D:\\Project\\宏基因组上游处理"
)
project <- "PRJNA834801"

taxonomy <- "species"

abundance <-
  read.xlsx(
    paste0(project, "_summary.xlsx"),
    sheet = taxonomy,
    rowNames = T
  )

abundance = floor(abundance)

##### metaData 读取
# metadata
metaData <-
  read.xlsx(
    paste0(project, "_summary.xlsx"),
    sheet = "metaData",
    rowNames = T
  )


# 微生物属级别丰度数据预处理
abundance <- abundance[, rownames(metaData)]

abundance_norm <-
  scale(abundance, center = FALSE, scale = colSums(abundance))

rownames(abundance_norm) <- str_split(rownames(abundance_norm),
                                      paste0(strsplit("species", "",)[[1]][1], "__"),
                                      simplify = T)[, 2]


data_species <- as.data.frame(t(abundance))

# 计算 Bray-Crutis距离
bray <- vegdist(data_species, method = "bray")

# 主坐标分析
pcoa = pcoa(bray)

# 取前两个主成分
pcoa_point = data.frame(pcoa$vectors[, c(1, 2)])

# 分组信息
pcoa_point$Case_status = metaData$Case_status
pcoa_point$Case_status = factor(pcoa_point$Case_status,
                               levels = unique(pcoa_point$Case_status))

# PERMANOVA 检验
set.seed(1)
data_species_div_dis <- adonis2(
  # ~ 分组
  data_species ~ metaData$Case_status,
  data = data_species,
  permutations = 999,
  method = "bray"
)
r2 <- round(data_species_div_dis$R2[1], 3)
f_values <- round(data_species_div_dis$F[1], 3)
p_value <- round(data_species_div_dis$`Pr(>F)`[1], 3)


# 坐标标签
xylab = paste(c("Axis.1 (", "Axis.2 ("),
              round(pcoa$values$Relative_eig[1:2] * 100, 2),
              "%)",
              sep = "")
linetype = 5
# 分组颜色设置,记得改成unique(分组)的数量
color <- pal_npg("nrc")(2)

# 绘制beta多样性聚类图
beta_diversity <-
  ggplot(data.frame(pcoa_point),
         # fill = 分组
         aes(x = Axis.1, y = Axis.2, fill = Case_status)) +
  # 置信区间绘制
  stat_ellipse(
    # color = 分组
    aes(x = Axis.1, y = Axis.2, color = Case_status),
    level = 0.95,
    show.legend = F,
    geom = "polygon",
    linetype = linetype,
    size = 1,
    alpha = 0.1
  ) +
  # 点绘制
  # color = 分组
  geom_point(aes(color = Case_status),
             size = 1) +
  # scale_shape_manual(values = c(21, 22, 23, 24, 25)) +
  theme_bw() +
  annotate(
    size = 0.1,
    geom = 'segment',
    y = Inf,
    yend = Inf,
    x = -Inf,
    xend = Inf
  ) +
  labs(x = xylab[1], y = xylab[2]) +
  scale_fill_manual(name = "Case_status", values = color) +
  scale_color_manual(name = "Case_status", values = color) +
  scale_x_continuous(limits = c(-1.0, 1.0)) +
  scale_y_continuous(limits = c(-.8, .8)) +
  theme(
    panel.border = element_rect(
      fill = NA ,
      size = 0.3,
      linetype = "solid"
    ),
    aspect.ratio = 1,
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # axis.line = element_line(color="gray", size = 0.1),
    text = element_text(family = "serif"),
    legend.title = element_text(face = 'bold', size = 12),
    legend.text = element_text(size = 8),
    axis.title = element_text(face = "bold", size = 12),
    axis.text = element_text(size = 8),
    plot.title = element_text(face = "bold", size = 10, hjust = 0.5)
  ) +
  # 附加文字注释
  annotate(
    "text",
    x = 0.6,
    y = -0.65,
    label = bquote("PERMANOVA"),
    parse = F,
    hjust = 0,
    size = 3,
    family = "serif"
  ) +
  annotate(
    "text",
    x = 0.6,
    y = -0.70,
    label = "R^2",
    parse = T,
    hjust = 0,
    size = 3,
    family = "serif"
  ) +
  annotate(
    "text",
    x = 0.67,
    y = -0.70,
    label = paste0(": ", r2),
    parse = F,
    hjust = 0,
    size = 3,
    family = "serif"
  ) +
  annotate(
    "text",
    x = 0.6,
    y = -0.75,
    label = "F value: ",
    hjust = 0,
    size = 3,
    family = "serif"
  ) +
  annotate(
    "text",
    x = 0.81,
    y = -0.75,
    label = f_values,
    hjust = 0,
    size = 3,
    family = "serif"
  ) +
  annotate(
    "text",
    x = 0.6,
    y = -0.80,
    label = "P value: ",
    hjust = 0,
    size = 3,
    family = "serif"
  ) +
  annotate(
    "text",
    x = 0.81,
    y = -0.80,
    label = p_value,
    hjust = 0,
    size = 3,
    color = "red",
    family = "serif"
  )

beta_diversity


ggsave(
  "β-Diversity.pdf",
  beta_diversity,
  width = 10,
  height = 8,
)





