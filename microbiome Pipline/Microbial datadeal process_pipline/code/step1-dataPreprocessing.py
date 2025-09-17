# -*- coding: utf-8 -*-
"""
Project         : 2022-01-10-1 01_癌症微生物数据库
File            : step1-dataPreprocessing.py
Creator         : Ping
Create time     : 2022/9/20 15:14
Software        : PyCharm
Introduction    :
"""

import os

import pandas as pd
import xlwings as xw

os.chdir(r"D:\Graduation_Design\1material\code")

# bracken 转换为 mpa 并合并
project = "GSE191328"
bracken_report_path = "raw_data/%s/brackenReport" % project

for file in [x for x in os.listdir(bracken_report_path) if "species" in x]:
    os.remove("%s/%s" % (bracken_report_path, file))

mpa_data_path = "mpa_data"
mpa_report_path = "mpa_data/mpa_report"
mpa_file_path = "%s/%s" % (mpa_report_path, project)

for path in [mpa_data_path, mpa_report_path, mpa_file_path]:
    if not os.path.exists(path):
        os.mkdir(path)


script1 = "mytools/brackenConvertTools/step1_kreport2mpa_WPP.py"

script2 = "mytools/brackenConvertTools/step2_combine_mpa_WPP.py"

command1 = "python %s -rp %s -op %s" % (script1,
                                        bracken_report_path, mpa_file_path)

os.system(command1)

command2 = "python %s -ip %s -o %s/%s_combinedMpa.txt" % (script2, mpa_file_path,
                                                          mpa_data_path, project)
os.system(command2)


# Data filtering, coexist as Excel


# Filter low-quality readings
# , mean_count=4, median_count=4):
def low_count_filter(data, min_count=4, min_percentage=0.1):
    taxonomy_names = data.index
    chosen_list = []
    for taxo in taxonomy_names:
        taxo_data = data.loc[taxo].tolist()
        if len([x for x in taxo_data if x >= min_count]) / len(taxo_data) >= min_percentage:
            chosen_list.append(taxo)
    return data.loc[chosen_list].copy()


# Obtain taxation level data
def get_taxonomy(data, taxo="s"):
    all_taxo = data.index.tolist()
    selected_taxo = [x for x in all_taxo if "%s__" % taxo in x.split("|")[-1]]
    selected_data = data.loc[selected_taxo].copy()
    # selected_data.index = [x.split("%s__" % taxo)[-1] for x in selected_data.index]
    selected_data_sorted = selected_data.loc[selected_data.mean(
        1).sort_values(ascending=False).index.tolist()]
    return selected_data_sorted


xlName = "analysis_data/%s_summary.xlsx" % project

if not os.path.exists(xlName):
    xlapp = xw.App(add_book=False)
    xlFile = xlapp.books.add()
    xlFile.save(xlName)
else:
    xlapp = xw.App(add_book=False)
    xlFile = xlapp.books.open(xlName)
# Create raw filtered


mpaData = pd.read_csv("mpa_data/%s_combinedMpa.txt" % project,
                      sep="\t", index_col=0)
mpaData_filtered = low_count_filter(mpaData)

all_sheets = xlFile.sheets
all_sheets_names = [x.name for x in all_sheets]
if "metaData" not in all_sheets_names:
    sheet = xlFile.sheets["Sheet1"]
    sheet.name = "metaData"
if "raw" not in all_sheets_names:
    xlFile.sheets.add("raw", after="metaData")
if "filtered" not in all_sheets_names:
    xlFile.sheets.add("filtered", after="raw")

# write metadata
metaData = pd.read_csv("analysis_data/%s_SraRunTable.txt" %
                       project, index_col=0)
sheet = xlFile.sheets["metaData"]
sheet.clear_contents()
sheet.range("A1").value = metaData
sheet.autofit()

sheet = xlFile.sheets["raw"]
sheet.clear_contents()
sheet.range("A1").value = mpaData
sheet.autofit()

sheet = xlFile.sheets["filtered"]
sheet.clear_contents()
sheet.range("A1").value = mpaData_filtered
sheet.autofit()

all_taxonomy = ["species", "genus", "family", "order", "class", "phylum"]
for i in range(len(all_taxonomy)):
    taxonomy = all_taxonomy[i]
    after = "filtered" if i == 0 else all_taxonomy[i - 1]
    if taxonomy not in all_sheets_names:
        xlFile.sheets.add(taxonomy, after=after)
    mpaData_selected = get_taxonomy(mpaData_filtered, taxonomy[0])
    sheet = xlFile.sheets[taxonomy]
    sheet.clear_contents()
    sheet.range("A1").value = mpaData_selected
    sheet.autofit()

sheet = xlFile.sheets["metaData"]
sheet.activate()
xlFile.save()
xlapp.quit()
