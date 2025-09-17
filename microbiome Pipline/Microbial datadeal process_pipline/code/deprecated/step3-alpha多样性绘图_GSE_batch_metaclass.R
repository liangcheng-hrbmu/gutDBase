library(pacman)
pacman::p_load(
  vegan, reshape2, dplyr, stringr, openxlsx
)

plot_alpha_diversity_data = function(
    accession,
    group_info,
    taxonomy,
    work_dir,
    host,
    meta_file) {
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

  if (length(meta_data$Sample.Name) != ncol(abundance)) {
    cat("跳过：", accession, "-", group_info, "下载数据不全\n")
    return(NULL)
  }

  abundance_norm = scale(abundance, center = FALSE, scale = colSums(abundance))
  rownames(abundance_norm) = str_split(
    rownames(abundance_norm),
    paste0(strsplit(taxonomy, "")[[1]][1], "__"),
    simplify = TRUE
  )[, 2]

  group_list = split(sub$sample, sub$value)
  sample_included = unlist(group_list)

  abundance = abundance[, sample_included, drop = FALSE]
  abundance_norm = abundance_norm[, sample_included, drop = FALSE]

  group_factor = rep(NA_character_, length(sample_included))
  names(group_factor) = sample_included
  for (grp in names(group_list)) {
    group_factor[group_list[[grp]]] = grp
  }

  group_factor = group_factor[sample_included]

  chao_ace = estimateR(t(abundance))
  chao1 = chao_ace[2, sample_included]
  shannon = diversity(t(abundance_norm), index = "shannon")[sample_included]
  simpson = diversity(t(abundance_norm), index = "simpson")[sample_included]

  data_plot = data.frame(
    sample = sample_included,
    chao1 = chao1,
    shannon = shannon,
    simpson = simpson,
    metaclass = group_factor,
    accession = accession
  )

  data_plot$id = seq_len(nrow(data_plot)) - 1
  data_plot = data_plot[, c("id", "sample", "chao1", "shannon", "simpson", "metaclass", "accession")]

  return(data_plot)
}

# 设置路径
meta_file = "D:/Project/gutDBase/metadata/hum_pie.csv"
work_dir = "D:/Project/gutDBase/bracken_summary_GSE"
host = "human"
taxonomy = "species"
output_dir = file.path("D:/Project/gutDBase/Alpha_Diversity_csv", host)
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

pie_meta = read.csv(meta_file, stringsAsFactors = FALSE)
accession_unique = unique(pie_meta[, "accession"])
combo_table = data.frame(accession = accession_unique, info = rep("metaclass", length(accession_unique)))

all_alpha_data = list()

for (i in seq_len(nrow(combo_table))) {
  acc = combo_table$accession[i]
  info = combo_table$info[i]
  cat("处理：", acc, "-", info, "\n")

  try({
    res = plot_alpha_diversity_data(
      accession = acc,
      group_info = info,
      taxonomy = taxonomy,
      work_dir = work_dir,
      host = host,
      meta_file = meta_file
    )
    if (!is.null(res)) {
      all_alpha_data[[length(all_alpha_data) + 1]] = res
    }
  })
}

if (length(all_alpha_data) > 0) {
  combined_alpha = do.call(rbind, all_alpha_data)
  write.csv(combined_alpha, file = file.path(output_dir, "combined_alpha_diversity.csv"), row.names = FALSE)
  cat("✅ 多样性指标数据已成功保存\n")
} else {
  cat("⚠️ 没有生成多样性数据\n")
}
