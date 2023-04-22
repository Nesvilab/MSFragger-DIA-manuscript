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


def translate_spectronaut_peptide(aa):
    aa.index = aa.index.str. \
        replace("[Acetyl (Protein N-term)]", "(UniMod:1)", regex=False).str. \
        replace("[Oxidation (M)]", "(UniMod:35)", regex=False).str. \
        replace("[Carbamidomethyl (C)]", "(UniMod:4)", regex=False).str. \
        replace("[Phospho (STY)]", "(UniMod:21)", regex=False).str. \
        replace("_", "", regex=False).str. \
        replace(".", "", regex=False)
    return aa


def translate_encyclopedia_peptides(aa):
    aa.index = aa.index.str. \
        replace("[+42.010565]", "(UniMod:1)", regex=False).str. \
        replace("[+15.994915]", "(UniMod:35)", regex=False).str. \
        replace("[+57.021464]", "(UniMod:4)", regex=False).str. \
        replace("[+79.966331]", "(UniMod:21)", regex=False)
    return aa


def process(input_df, overlap_df, other_union_df, tool_name):
    input_overlap = input_df.loc[overlap_df]
    tt = input_overlap.dropna(axis=0, thresh=3)
    input_overlap_cv = pd.DataFrame({"cv": np.nanstd(tt, 1) * 100 / np.nanmean(tt, 1), "tool": tool_name, "type": "overlap"})

    input_unique = input_df.loc[~input_df.index.isin(other_union_df)]
    tt = input_unique.dropna(axis=0, thresh=3)
    input_unique_cv = pd.DataFrame({"cv": np.nanstd(tt, 1) * 100 / np.nanmean(tt, 1), "tool": tool_name, "type": "unique"})

    return input_overlap_cv, input_unique_cv


os.chdir(r"G:\Dropbox\papers_Fengchao\msfragger_dia\script\results\2018-HeLa")

# plot a
spectronaut_path = r"spectronaut\modified_sequence_maxlfq.tsv"
encyclopedia_path = r"encyclopedia\encyclopedia_quant.elib.peptides.txt"
diann_path = r"diann\modified_sequence_maxlfq.tsv"
maxdia_path = r"maxdia\combined\txt\modificationSpecificPeptides.txt"
msfraggerdia_path = r"msfraggerdia\modified_sequence_maxlfq.tsv"

spectronaut = pd.read_csv(spectronaut_path, sep="\t", index_col=0, na_values="NA", header=0, usecols=["Unnamed: 0", "23aug2017_hela_serum_timecourse_pool_wide_001", "23aug2017_hela_serum_timecourse_pool_wide_002", "23aug2017_hela_serum_timecourse_pool_wide_003"])

encyclopedia = pd.read_csv(encyclopedia_path, sep="\t", index_col=0, na_values=0, header=0, usecols=["Peptide", "23aug2017_hela_serum_timecourse_pool_wide_001.mzML", "23aug2017_hela_serum_timecourse_pool_wide_002.mzML", "23aug2017_hela_serum_timecourse_pool_wide_003.mzML"])

diann = pd.read_csv(diann_path, sep="\t", index_col=0, na_values="NA", header=0, usecols=["Unnamed: 0", "23aug2017_hela_serum_timecourse_pool_wide_001", "23aug2017_hela_serum_timecourse_pool_wide_002", "23aug2017_hela_serum_timecourse_pool_wide_003"])

maxdia = pd.read_csv(maxdia_path, sep="\t", index_col=0, na_values=[0, "", "NaN"], header=0, usecols=["Sequence", "Intensity 23aug2017_hela_serum_timecourse_pool_wide_001", "Intensity 23aug2017_hela_serum_timecourse_pool_wide_002", "Intensity 23aug2017_hela_serum_timecourse_pool_wide_003", "Reverse", "Potential contaminant"])

msfraggerdia = pd.read_csv(msfraggerdia_path, sep="\t", index_col=0, na_values="NA", header=0, usecols=["Unnamed: 0", "23aug2017_hela_serum_timecourse_pool_wide_001", "23aug2017_hela_serum_timecourse_pool_wide_002", "23aug2017_hela_serum_timecourse_pool_wide_003"])

spectronaut.dropna(how="all", inplace=True)
encyclopedia.dropna(how="all", inplace=True)
diann.dropna(how="all", inplace=True)

maxdia = maxdia.loc[pd.isna(maxdia["Reverse"]) & pd.isna(maxdia["Potential contaminant"])]
maxdia.drop(labels=["Reverse", "Potential contaminant"], axis=1, inplace=True)
maxdia.dropna(how="all", inplace=True)

msfraggerdia.dropna(how="all", inplace=True)

spectronaut = translate_spectronaut_peptide(spectronaut)
encyclopedia = translate_encyclopedia_peptides(encyclopedia)

spectronaut_count = pd.DataFrame({"count": spectronaut.count(axis=0), "tool": "Spectronaut 17"})
publication_count = pd.DataFrame({"count": [65219, 65675, 65845], "tool": "Searle et al. 2020"})
encyclopedia_count = pd.DataFrame({"count": encyclopedia.count(axis=0), "tool": "EncyclopeDIA"})
diann_count = pd.DataFrame({"count": diann.count(axis=0), "tool": "DIA-NN\nlib-free"})
maxdia_count = pd.DataFrame({"count": maxdia.count(axis=0), "tool": "MaxDIA"})
msfraggerdia_count = pd.DataFrame({"count": msfraggerdia.count(axis=0), "tool": "FP-MSF"})

average_counts = pd.DataFrame({"Spectronaut 17": spectronaut_count["count"].mean(), "Searle et al. 2020": publication_count["count"].mean(), "EncyclopeDIA": encyclopedia_count["count"].mean(), "DIA-NN\nlib-free": diann_count["count"].mean(), "FP-MSF": msfraggerdia_count["count"].mean()}, index=[0])

run_counts = pd.concat([spectronaut_count, publication_count, encyclopedia_count, diann_count, msfraggerdia_count])

sns_plot = sns.barplot(data=average_counts, palette=[palette_dict["spectronaut"], palette_dict["publication"], palette_dict["encyclopedia"], palette_dict["diann"], palette_dict["msfraggerdia"]])
sns_plot = sns.swarmplot(x="tool", y="count", data=run_counts, color="white", edgecolor="black", size=10, linewidth=2)
sns_plot.set(xlabel=None, ylabel="quantified peptides")
sns_plot.set_xticklabels(sns_plot.get_xticklabels(), fontsize=9)
sns_plot.set_xticklabels(sns_plot.get_xticklabels())
sns_plot.figure.savefig("2018-HeLa_a.pdf")
plt.figure()


# plot a S
average_counts = pd.DataFrame({"Spectronaut 17": spectronaut_count["count"].mean(), "Searle et al. 2020": publication_count["count"].mean(), "EncyclopeDIA": encyclopedia_count["count"].mean(), "DIA-NN\nlib-free": diann_count["count"].mean(), "MaxDIA": maxdia_count["count"].mean(), "FP-MSF": msfraggerdia_count["count"].mean()}, index=[0])

run_counts = pd.concat([spectronaut_count, publication_count, encyclopedia_count, diann_count, maxdia_count, msfraggerdia_count])

sns_plot = sns.barplot(data=average_counts, palette=[palette_dict["spectronaut"], palette_dict["publication"], palette_dict["encyclopedia"], palette_dict["diann"], palette_dict["maxdia"], palette_dict["msfraggerdia"]])
sns_plot = sns.swarmplot(x="tool", y="count", data=run_counts, color="white", edgecolor="black", size=10, linewidth=2)
sns_plot.set(xlabel=None, ylabel="quantified peptides")
sns_plot.set_xticklabels(sns_plot.get_xticklabels(), fontsize=8)
sns_plot.figure.savefig("2018-HeLa_a_S.pdf")
plt.figure()


# plot b
all_overlap = spectronaut.index.intersection(encyclopedia.index).intersection(diann.index).intersection(msfraggerdia.index)

other_union = encyclopedia.index.union(diann.index).union(msfraggerdia.index)
spectronaut_overlap_cv, spectronaut_unique_cv = process(spectronaut, all_overlap, other_union, "Spectronaut 17")

other_union = spectronaut.index.union(diann.index).union(msfraggerdia.index)
encyclopedia_overlap_cv, encyclopedia_unique_cv = process(encyclopedia, all_overlap, other_union, "EncyclopeDIA")

other_union = spectronaut.index.union(encyclopedia.index).union(msfraggerdia.index)
diann_overlap_cv, diann_unique_cv = process(diann, all_overlap, other_union, "DIA-NN\nlib-free")

other_union = spectronaut.index.union(encyclopedia.index).union(diann.index)
msfraggerdia_overlap_cv, msfraggerdia_unique_cv = process(msfraggerdia, all_overlap, other_union, "FP-MSF")

x = pd.concat([spectronaut_overlap_cv, encyclopedia_overlap_cv, diann_overlap_cv, msfraggerdia_overlap_cv, spectronaut_unique_cv, encyclopedia_unique_cv, diann_unique_cv, msfraggerdia_unique_cv])

sns_plot = sns.boxplot(x="tool", y="cv", hue="type", hue_order=["overlap", "unique"], data=x, linewidth=1, showfliers=False)
sns_plot.set(xlabel=None, ylabel="coefficient of variation (%)")
plt.ylim(0, 100)
sns_plot.figure.savefig("2018-HeLa_b.pdf")
plt.figure()
