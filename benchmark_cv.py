#  Licensed to the Apache Software Foundation (ASF) under one
#  or more contributor license agreements.  See the NOTICE file
#  distributed with this work for additional information
#  regarding copyright ownership.  The ASF licenses this file
#  to you under the Apache License, Version 2.0 (the
#  "License"); you may not use this file except in compliance
#  with the License.  You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
#  Unless required by applicable law or agreed to in writing,
#  software distributed under the License is distributed on an
#  "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
#  KIND, either express or implied.  See the License for the
#  specific language governing permissions and limitations
#  under the License.
import os

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

sns.set(font="Arial")
sns.set_theme(style="ticks")

palette_dict = {
    "publication": sns.color_palette()[9],
    "diaumpire": sns.color_palette()[6],
    "diaumpire hybrid": sns.color_palette()[8],
    "spectronaut": sns.color_palette()[4],
    "encyclopedia": sns.color_palette()[5],
    "diann": sns.color_palette()[3],
    "maxdia": sns.color_palette()[2],
    "msfraggerdia": sns.color_palette()[0],
    "msfragger hybrid": sns.color_palette()[1]
}


def pick_normalize_ecoli_precursors(aa, annotation2_path, is_spectronaut):
    if is_spectronaut:
        annotation2 = pd.read_csv(annotation2_path, sep="\t", index_col=None, na_values="NA", header=0, usecols=["PG.ProteinGroups", "PEP.StrippedSequence", "EG.ModifiedSequence", "FG.Charge"])
        annotation2["Precursor.Id"] = annotation2["EG.ModifiedSequence"] + annotation2["FG.Charge"].apply(str)
        annotation2.rename(columns={"PG.ProteinGroups": "Protein.Group", "PEP.StrippedSequence": "Stripped.Sequence"}, inplace=True)
    else:
        annotation2 = pd.read_csv(annotation2_path, sep="\t", index_col=None, na_values="NA", header=0, usecols=["Precursor.Id", "Protein.Group", "Stripped.Sequence"])
    annotation2.drop_duplicates(subset="Precursor.Id", keep="first", inplace=True)
    annotation2 = annotation2.loc[:, ["Precursor.Id", "Protein.Group", "Stripped.Sequence"]]
    annotation2["Protein.Group"] = annotation2["Protein.Group"].str.split(";", expand=True, regex=False)[0]
    annotation2 = annotation2.merge(annotation, how="left", left_on="Protein.Group", right_on="Entry")
    annotation2.set_index("Precursor.Id", drop=True, inplace=True)
    aa.index.name = "Precursor.Id"
    aa = aa.merge(annotation2, how="left", on="Precursor.Id")
    aa = aa[aa["Organism"] == "Escherichia coli (strain K12)"]
    aa = aa[~aa["Stripped.Sequence"].isin(human_ecoli_overlap_peptide)]
    aa.drop(["Protein.Group", "Stripped.Sequence", "Organism"], axis=1, inplace=True)
    median_column = aa.median(axis=0, skipna=True)
    median_all = median_column.median(axis=0, skipna=True)
    aa = aa * median_all / median_column
    return aa


def translate_spectronaut_peptides(aa):
    aa.index = aa.index.str.\
        replace("[Acetyl (Protein N-term)]", "(UniMod:1)", regex=False).str.\
        replace("[Oxidation (M)]", "(UniMod:35)", regex=False).str.\
        replace("[Carbamidomethyl (C)]", "(UniMod:4)", regex=False).str.\
        replace("[Phospho (STY)]", "(UniMod:21)", regex=False).str.\
        replace("_", "", regex=False).str. \
        replace(".", "", regex=False)
    return aa


def process(input_df, tool_name):
    xx = input_df.dropna(axis=0, thresh=4)
    input_cv = pd.DataFrame({"cv": np.nanstd(xx, 1) * 100 / np.nanmean(xx, 1), "tool": tool_name})
    return input_cv


os.chdir(r"G:\Dropbox\papers_Fengchao\msfragger_dia\script\results\benchmark")

annotation = pd.read_csv(r"G:\Dropbox\papers_Fengchao\msfragger_dia\script\code\MSFragger-DIA-manuscript\uniprot.tab", sep="\t", index_col=0, na_values="NA", header=0, usecols=["Entry", "Organism"])

human_ecoli_overlap_peptide = []
with open("human_ecoli_overlap_peptides.txt") as f:
    for line in f.readlines():
        human_ecoli_overlap_peptide.append(line.strip())

# precursor level
spectronaut_14_path = r"spectronaut\14\\"
spectronaut_17_path = r"spectronaut\17\\"
diann_path = r"diann\\"
msfraggerdia_path = r"msfraggerdia\\"
msfraggerdia_hybrid_path = r"msfraggerdiadda\\"

spectronaut_14 = pd.read_csv(spectronaut_14_path + "precursor_maxlfq.tsv", sep="\t", index_col=0, na_values="NA", header=0)
spectronaut_17 = pd.read_csv(spectronaut_17_path + "precursor_maxlfq.tsv", sep="\t", index_col=0, na_values="NA", header=0)
diann = pd.read_csv(diann_path + "precursor_maxlfq.tsv", sep="\t", index_col=0, na_values="NA", header=0)
msfraggerdia = pd.read_csv(msfraggerdia_path + "precursor_maxlfq.tsv", sep="\t", index_col=0, na_values="NA", header=0)
msfraggerdia_hybrid = pd.read_csv(msfraggerdia_hybrid_path + "precursor_maxlfq.tsv", sep="\t", index_col=0, na_values="NA", header=0)

spectronaut_14.dropna(how="all", inplace=True)
spectronaut_17.dropna(how="all", inplace=True)
diann.dropna(how="all", inplace=True)
msfraggerdia.dropna(how="all", inplace=True)
msfraggerdia_hybrid.dropna(how="all", inplace=True)

spectronaut_14 = pick_normalize_ecoli_precursors(spectronaut_14, spectronaut_14_path + "20230227_181426_directDIA_LymphEcoli_Report.tsv", True)
spectronaut_17 = pick_normalize_ecoli_precursors(spectronaut_17, spectronaut_17_path + "20230227_180032_SN17_RealDilutionSeries_Report.tsv", True)
diann = pick_normalize_ecoli_precursors(diann, diann_path + "diann-output.tsv", False)
msfraggerdia = pick_normalize_ecoli_precursors(msfraggerdia, msfraggerdia_path + "diann-output.tsv", False)
msfraggerdia_hybrid = pick_normalize_ecoli_precursors(msfraggerdia_hybrid, msfraggerdia_hybrid_path + "diann-output.tsv", False)

spectronaut_14 = translate_spectronaut_peptides(spectronaut_14)
spectronaut_17 = translate_spectronaut_peptides(spectronaut_17)

spectronaut_14.drop_duplicates(keep="first", inplace=True)
spectronaut_17.drop_duplicates(keep="first", inplace=True)
diann.drop_duplicates(keep="first", inplace=True)
msfraggerdia.drop_duplicates(keep="first", inplace=True)
msfraggerdia_hybrid.drop_duplicates(keep="first", inplace=True)

conditions = ["1-06"]
condition_patterns = [r"Lymph_Ecoli_1-6_\d+"]

for idx, p in enumerate(condition_patterns):
    spectronaut_14_sub = spectronaut_14.filter(regex=p, axis=1)
    spectronaut_17_sub = spectronaut_17.filter(regex=p, axis=1)
    diann_sub = diann.filter(regex=p, axis=1)
    msfraggerdia_sub = msfraggerdia.filter(regex=p, axis=1)
    msfraggerdia_hybrid_sub = msfraggerdia_hybrid.filter(regex=p, axis=1)

    spectronaut_14_cv = process(spectronaut_14_sub, "Spectronaut\n14")
    spectronaut_17_cv = process(spectronaut_17_sub, "Spectronaut\n17")
    diann_cv = process(diann_sub, "DIA-NN\nlib-free")
    msfraggerdia_cv = process(msfraggerdia_sub, "FP-MSF")
    msfraggerdia_hybrid_cv = process(msfraggerdia_hybrid_sub, "FP-MSF hybrid")

    x = pd.concat([spectronaut_14_cv, spectronaut_17_cv, diann_cv, msfraggerdia_cv, msfraggerdia_hybrid_cv])

    plt.figure(figsize=(6, 6), dpi=300)
    sns_plot = sns.violinplot(x="tool", y="cv", data=x, inner=None, linewidth=0.5, width=1, cut=0, palette=[palette_dict["spectronaut"], palette_dict["spectronaut"], palette_dict["diann"], palette_dict["maxdia"], palette_dict["msfraggerdia"]])
    sns.boxplot(x="tool", y="cv", data=x, showfliers=False, width=0.3, linewidth=1, boxprops={'zorder': 2, 'facecolor': 'white'}, ax=sns_plot)
    sns_plot.set_ylim([0, 100])
    sns_plot.set(xlabel=None, ylabel="coefficient of variation (%)")
    sns_plot.figure.savefig("benchmark_cv_precursor_" + conditions[idx] + ".pdf")
