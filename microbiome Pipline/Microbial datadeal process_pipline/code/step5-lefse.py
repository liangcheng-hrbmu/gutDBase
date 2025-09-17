from collections import Counter
import pandas as pd
import os
import seaborn as sns

os.chdir(r"D:\Graduation_Design\1material\code")
project = "GSE191328"

abundance = pd.read_excel("analysis_data/%s_summary.xlsx" %
                          project, sheet_name="filtered", index_col=0)

metaData = pd.read_excel("analysis_data/%s_summary.xlsx" %
                         project, sheet_name="metaData", index_col=0).T

metaData = metaData[metaData.columns[metaData.loc["Time_point"].isin([
                                                                     "baseline"])]]

abundance = abundance[metaData.columns]

lefse_input = pd.concat(
    [metaData.columns.to_frame().T, metaData.loc[["diagnosis"]], abundance])

lefse_input.index.name = "SampleID"

lefse_input.to_csv("lefse/%s_lefse_input_diagnosis.txt" %
                   project, sep="\t", header=False)

# Server

lefse_step1 = "lefse_format_input.py %s_lefse_input_diagnosis.txt %s_lefse_diagnosis.in -c 2 -u -1 -o 1000000" % (
    project, project)

lefse_step2 = "lefse_run.py %s_lefse_diagnosis.in %s_lefse_diagnosis.res -l 1.5" % (
    project, project)

lefse_step3 = "lefse_plot_cladogram.py %s_lefse_diagnosis.res %s_lefse_cladogram_diagnosis.pdf --format pdf --abrv_stop_lev 7" % (
    project, project)


# Preprocess the results and draw a circular diagram of the evolutionary tree

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

taxo_dict = dict(zip(['p', 'c', 'o', 'f', 'g', 's'],
                     ['phylum', 'class', 'order', 'family', 'genus', 'species']))

# Before reading, it is necessary to save% s_lfse_input-diagnosis.xls as% s_lfse_input-diagnosis.xlsx
lefse_result = pd.read_excel("lefse/%s_lefse_input_diagnosis.xlsx" %
                             project, sheet_name="%s_lefse_input_diagnosis" % project, header=None)

lefse_result.columns = [
    "Taxonomy", "Log(max(meanAbundance))", "group", "LDA score", "P_value"]

lefse_result["Taxonomy"] = [x.replace("s___Mycobacterium__stephanolepidis", "s__Mycobacterium_stephanolepidis") for x in
                            lefse_result["Taxonomy"]]


lefse_result = lefse_result.dropna()

lefse_result.insert(0, column="taxon",
                    value=[x.split(".")[-1][3:] for x in lefse_result["Taxonomy"].tolist()])

lefse_result.insert(0, column="taxonomy_level",
                    value=[taxo_dict[x.split(".")[-1][0]] for x in lefse_result["Taxonomy"].tolist()])

lefse_input.index = [multiple_replace(
    x, replace_dict) for x in lefse_input.index]


lefse_input.loc[[x for x in lefse_input.index if "s__" in x]].to_csv("lefse/%s_heatmap_data_diagnosis.txt" % project,
                                                                     sep="\t")


abundance_clad = abundance.copy()

abundance_clad.index = [multiple_replace(
    x, replace_dict) for x in abundance_clad.index]

idmap = dict(zip([x.split("s__")[-1]
             for x in abundance_clad.index], abundance_clad.index))

otuTab_data = lefse_input.loc[lefse_input.index[:1].tolist() +
                              [idmap[x] for x in lefse_result.loc[lefse_result["taxonomy_level"] == "species"][
                                  "taxon"].tolist()]]


keys = ['p', 'c', 'o', 'f', 'g', 's']
value = []
for tax in [idmap[x] for x in lefse_result.loc[lefse_result["taxonomy_level"] == "species"]["taxon"].tolist()]:
    tax_dict = dict(zip([x[0] for x in tax.split("|")],
                    [x[3:] for x in tax.split("|")]))
    value.append([tax_dict[x] if x in tax_dict else "__None" for x in keys])

taxonomy_data = pd.DataFrame(
    columns=['phylum', 'class', 'order', 'family', 'genus', 'species'], data=value)

taxonomy_data.index = taxonomy_data["species"]
taxonomy_data.index.name = "biomarkers"
taxonomy_data["tree_bone"] = [".".join(x) for x in
                              taxonomy_data[['phylum', 'class', 'order', 'family', 'genus', 'species']].values]

taxonomy_data["group"] = lefse_result.loc[lefse_result["taxonomy_level"]
                                          == "species"]["group"].tolist()

taxonomy_data["LDA score"] = lefse_result.loc[lefse_result["taxonomy_level"]
                                              == "species"]["LDA score"].tolist()

taxonomy_data.to_csv(
    "lefse/lefse_result_taxonomy_data_%s_diagnosis.txt" % project, sep="\t")

taxonomy_data.groupby(["group", "phylum"]).count()

heatmap_data = abundance_clad.loc[[idmap[x] for x in taxonomy_data["species"]]]
heatmap_data.to_csv("lefse/heatmap_abundanceData_diagnosis.txt", sep="\t")


Counter(taxonomy_data["group"])
Counter(lefse_input.loc["diagnosis"])

with open("lefse/tree_backbone_diagnosis.txt", "w") as fw:
    fw.write("\n".join(taxonomy_data["tree_bone"]))

otuTab = otuTab_data.T

otuTab[otuTab.columns[1:]] = otuTab[otuTab.columns[1:]].astype(float)

abundance_statistics = otuTab.groupby(otuTab.columns[0]).mean()

abundance_statistics_1 = abundance_statistics / abundance_statistics.sum()
abundance_statistics_1.columns = [
    x.split("s__")[-1] for x in abundance_statistics_1.columns]

with open("lefse/global_annotation.txt") as fr:
    global_annotation = fr.read()

phys_all = sorted(set(taxonomy_data["phylum"].tolist()))
phy_color_marker = dict(
    zip(phys_all, sns.hls_palette(len(phys_all), l=.5, s=0.7).as_hex()))
phy_color_background = dict(
    zip(phys_all, sns.hls_palette(len(phys_all), l=.7, s=0.5).as_hex()))
tax_phy_dict = dict(zip(taxonomy_data["tree_bone"], taxonomy_data["phylum"]))

with open("lefse/graphlan_annotate_diagnosis.txt", "w") as fw:
    global_annotation += "\n\n### phylum level annotation_background_color\n"
    # Background color annotation
    for phy in phys_all:
        global_annotation += "%s\tannotation_background_color\t%s\n" % (
            phy, phy_color_marker[phy])

    # Node annotation
    global_annotation += "\n\n###annotation\n"
    for tax in sorted(tax_phy_dict):
        global_annotation += "%s\tannotation\t*\n" % tax + \
                             "%s\tannotation_rotation\t90\n" % tax + \
                             "%s\tannotation_background_color\t%s\n" % (
                                 tax, phy_color_marker[tax_phy_dict[tax]])
    global_annotation += "\n\n###clade\n"
    for phy in phys_all:
        global_annotation += "%s*\tclade_marker_color\t%s\n" % (
            phy, phy_color_marker[phy])
        global_annotation += "%s\tclade_marker_size\t%s\n" % (phy, 100)
    global_annotation += "None\tclade_marker_size\t0\n"
    global_annotation += "\n\n###legend\n"
    for phy in phys_all:
        global_annotation += "%s\tclade_marker_color\t%s\n" % (
            "p_%s" % phy, phy_color_marker[phy])
        global_annotation += "%s\tclade_marker_size\t%s\n" % ("p_%s" % phy, 18)

    # Environmental annotation
    global_annotation += "\n\n###ring annotation 4\n"
    ring_group_color_dict = dict(
        zip(sorted(set(taxonomy_data["group"])), ["#D35800", "#008B44"]))
    for tax, group in zip(taxonomy_data["tree_bone"], taxonomy_data["group"]):
        global_annotation += "%s\tring_color\t4\t%s\n" % (
            tax, ring_group_color_dict[group])

    global_annotation += "\n\n###ring annotation 3\n"
    alpha_mitigant = 0.05
    ring_alphas = (taxonomy_data["LDA score"] - taxonomy_data["LDA score"].min() + alpha_mitigant) / (
        taxonomy_data["LDA score"].max() - taxonomy_data["LDA score"].min() + alpha_mitigant)
    ring_score_color_dict = dict(
        zip(sorted(set(taxonomy_data["group"])), ["#FF0000", "#2882D3"]))
    for tax, group, alpha in zip(taxonomy_data["tree_bone"], taxonomy_data["group"], ring_alphas):
        global_annotation += "%s\tring_color\t3\t%s\n" % (
            tax, ring_score_color_dict[group])
        global_annotation += "%s\tring_alpha\t3\t%s\n" % (tax, alpha)
    global_annotation += "\n\n###ring legend annotation\n"
    for group in ring_group_color_dict:
        global_annotation += "%s\tclade_marker_color\t%s\n" % (
            "group_%s" % group, ring_group_color_dict[group])
        global_annotation += "%s\tclade_marker_shape\t%s\n" % (
            "group_%s" % group, "s")
        global_annotation += "%s\tclade_marker_size\t%s\n" % (
            "group_%s" % group, 18)

    global_annotation += "\n\n###ring annotation 2 group1\n"
    group = sorted(set(taxonomy_data["group"]))[0]
    global_annotation += " ring_label\t2\t%s\n" % group
    for tax, species in zip(taxonomy_data["tree_bone"], taxonomy_data["species"]):
        ring_alpha = abundance_statistics_1[species][group]
        global_annotation += "%s\tring_color\t2\t%s\n" % (tax, "#C02A00")
        global_annotation += "%s\tring_alpha\t2\t%s\n" % (tax, ring_alpha)
        global_annotation += "%s\tring_shape\t2\t%s\n" % (tax, "v")

    global_annotation += "\n\n###ring annotation 1 group1\n"
    group = sorted(set(taxonomy_data["group"]))[1]
    global_annotation += " ring_label\t1\t%s\n" % group
    for tax, species in zip(taxonomy_data["tree_bone"], taxonomy_data["species"]):
        ring_alpha = abundance_statistics_1[species][group]
        global_annotation += "%s\tring_color\t1\t%s\n" % (tax, "#0054C0")
        global_annotation += "%s\tring_alpha\t1\t%s\n" % (tax, ring_alpha)
        global_annotation += "%s\tring_shape\t1\t%s\n" % (tax, "^")

    fw.write(global_annotation)


lefse_result_genus = lefse_result.loc[lefse_result['taxonomy_level'] == 'genus'][[
    'taxon', 'group', 'LDA score']]
lefse_result_genus.to_csv("genus.txt", sep='\t')
lefse_result_species = lefse_result.loc[lefse_result['taxonomy_level'] == 'species'][[
    'taxon', 'group', 'LDA score']]
lefse_result_species.to_csv("species.txt", sep='\t')
