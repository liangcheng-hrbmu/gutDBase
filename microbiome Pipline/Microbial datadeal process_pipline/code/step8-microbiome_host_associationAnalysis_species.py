# -*- coding: utf-8 -*- 
"""
Project         : 2021-10-09-6 01_Sars2_Microbiome_GSE169687
File            : step8-microbiome_host_associationAnalysis.py
Creator         : Ping
Create time     : 2022/9/25 17:20
Software        : PyCharm
Introduction    :
"""

import os
from collections import Counter
import pandas as pd
from scipy.stats import pearsonr
from tqdm import tqdm

os.chdir(r"D:\Graduation_Design\1material\code")
project = "GSE191328"


def multiple_replace(string, *params):
    if len(params) == 1:
        if type(params[0]) is dict:
            replace_dict = params[0]
            for key in replace_dict:
                string = string.replace(key, replace_dict[key])
            return string
    if len(params) == 2:
        old = params[0]
        new = params[1]
        if type(old) is not list:
            if type(old) is str:
                old = [old]
            else:
                raise TypeError("data type error")
        if type(new) is not str:
            raise TypeError("data type error")
        for key in old:
            string = string.replace(key, new)
        return string


replace_dict = {
    "-": "_",
    ".": "_",
    "/": "_",
    "(": "_",
    ")": "_",
    "[": "",
    "]": "",
    "'": "",
}
# metaData = pd.read_csv("analysis_data/GSE169687_SraRunTable_orig.txt")
#
# gsm_srr_convert_dict = dict(zip(metaData["Sample Name"], metaData["Run"]))
#
# 合并宿主基因文件
# file_list = sorted(os.listdir("host_gene_data/GSE169687_RAW"))
# combined_data = pd.DataFrame()
# for file in file_list:
#     print("%3d of %3d %s" % (file_list.index(file) + 1, len(file_list), file))
#     sample_name = gsm_srr_convert_dict[file.split("_")[0]]
#
#     sample_data = pd.read_csv("host_gene_data/GSE169687_RAW/%s" % file,
#                               sep="\t",
#                               comment="#",
#                               index_col= 0)
#     combined_data = pd.concat([combined_data, sample_data[sample_data.columns[-1]]],
#                               axis=1)
#     combined_data.rename(columns={sample_data.columns[-1]: sample_name}, inplace=True)
#
# combined_data.index.name = "Geneid"
#
# # 转换ID
# convert_data = pd.read_csv("host_gene_data/ensembl_data.df.csv",
#                            sep="\t", dtype=str, index_col=0)
#
# # 去空行
# convert_data = convert_data.dropna(subset="external_gene_name")
#
# # 保留蛋白编码
# convert_data = convert_data.loc[convert_data["gene_biotype"] == "protein_coding"]
#
# # 取交集并转换
# combined_data_convert = combined_data.loc[sorted(set(combined_data.index.tolist()).
#                                                  intersection(set(convert_data["ensembl_gene_id"].tolist())))]
# convert_data = convert_data.loc[combined_data_convert.index]
#
#
# combined_data_convert.insert(0, column="gene_name", value= convert_data["external_gene_name"])
#
# combined_data_convert_1 = combined_data_convert.groupby("gene_name").mean()
#
# combined_data_convert_1.to_csv("host_gene_data/gene_expression_reads.txt", sep="\t")

abundance = pd.read_excel("analysis_data/%s_summary.xlsx" % project, sheet_name="species", index_col=0)

abundance.index = [multiple_replace(x.split("|s__")[-1], replace_dict) for x in abundance.index]

metaData = pd.read_excel("analysis_data/%s_summary.xlsx" % project, sheet_name="metaData", index_col=0)

metaData = metaData[metaData["Time_point"].isin(["baseline"])]

gene_exp = pd.read_csv("host_gene_data/GSE191328_exp_TPM_gene_name.txt", sep="\t", index_col=0)[
    metaData.index]

# gene_exp = gene_exp/gene_exp.sum()
# gene_exp.sum()

lefse_result = pd.read_csv("lefse/%s_lefse_input_diagnosis.xls" % project, sep="\t", header=None)

lefse_result.columns = ["Taxonomy", "Log(max(meanAbundance))", "group", "LDA score", "P_value"]
lefse_result["Taxonomy"] = [x.replace("s___Mycobacterium__stephanolepidis", "s__Mycobacterium_stephanolepidis") for x in
                            lefse_result["Taxonomy"]]
lefse_result.index = lefse_result["Taxonomy"]
lefse_result = lefse_result.dropna()
lefse_result = lefse_result[lefse_result["group"].isin(["UC", "CD"])]
lefse_result = lefse_result.loc[[x for x in lefse_result.index if ".s__" in x]]

biomarkers_treatment = [x.split(".s__")[-1] for x in lefse_result.index]

microbiome_abu = abundance.loc[biomarkers_treatment][metaData.index].copy()

# calculate correction score between gene_exp and microbiome_abu

gene_exp_norm = gene_exp / gene_exp.sum()
microbiome_abu_norm = microbiome_abu / microbiome_abu.sum()

cor_result = []
for spe in microbiome_abu.index:
    print(spe)
    for gene in gene_exp.index:
        corr = pearsonr(microbiome_abu_norm.loc[spe],
                        gene_exp_norm.loc[gene])
        cor_result.append([spe, gene, corr[0], corr[1]])

cor_result_data = pd.DataFrame(cor_result, columns=['species', 'gene', 'cor', 'p_values'])

cor_result_data = cor_result_data.dropna()


cor_result_data.to_csv("microbiome_host_association_result/cor_result_diagnosis.txt", sep="\t", index=False)

