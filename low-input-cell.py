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

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os

sns.set(font="Arial")
sns.set_theme(style="ticks")

palette_dict = {
    "publication": sns.color_palette()[9],
    "diaumpire": sns.color_palette()[6],
    "diaumpire hybrid": sns.color_palette()[8],
    "spectronaut": sns.color_palette()[4],
    "diann": sns.color_palette()[3],
    "maxdia": sns.color_palette()[2],
    "msfraggerdia": sns.color_palette()[0],
    "msfragger hybrid": sns.color_palette()[1]
}


def process(input_df, overlap_df, other_union_df, tool_name):
    input_overlap = input_df.loc[overlap_df]
    tt = input_overlap.dropna(axis=0, thresh=3)
    input_overlap_cv = pd.DataFrame({"cv": np.nanstd(tt, 1) * 100 / np.nanmean(tt, 1), "tool": tool_name, "type": "overlap"})

    input_unique = input_df.loc[~input_df.index.isin(other_union_df)]
    tt = input_unique.dropna(axis=0, thresh=3)
    input_unique_cv = pd.DataFrame({"cv": np.nanstd(tt, 1) * 100 / np.nanmean(tt, 1), "tool": tool_name, "type": "unique"})

    return input_overlap_cv, input_unique_cv


os.chdir(r"G:\Dropbox\papers_Fengchao\msfragger_dia\script\results")


# plot protein barplot
diaumpire_path_list = [r"low-input-cell_0.75ng_1.5ng\diaumpire\protein_maxlfq.tsv", r"low-input-cell_7.5ng_1ug\diaumpire\protein_maxlfq.tsv"]
diann_path_list = [r"low-input-cell_0.75ng_1.5ng\diann\protein_maxlfq.tsv", r"low-input-cell_7.5ng_1ug\diann\protein_maxlfq.tsv"]
maxdia_path_list = [r"low-input-cell_0.75ng_1.5ng\maxdia\combined\txt\proteinGroups.txt", r"low-input-cell_7.5ng_1ug\maxdia\combined\txt\proteinGroups.txt"]
msfraggerdia_path_list = [r"low-input-cell_0.75ng_1.5ng\msfraggerdia\protein_maxlfq.tsv", r"low-input-cell_7.5ng_1ug\msfraggerdia\protein_maxlfq.tsv"]

experiment_list = ["0.75 ng", "7.5 ng"]
msfraggerdia_cols = [["Unnamed: 0", "20210430_PC9-750pg_60K-IT118-DIA-IW10_01", "20210430_PC9-750pg_60K-IT118-DIA-IW10_02", "20210430_PC9-750pg_60K-IT118-DIA-IW10_03"],
                     ["Unnamed: 0", "20210430_PC9-7500pg_60K-IT118-DIA-IW10_01", "20210430_PC9-7500pg_60K-IT118-DIA-IW10_02", "20210430_PC9-7500pg_60K-IT118-DIA-IW10_03"]]

for expIdx in range(len(experiment_list)):
    diaumpire = pd.read_csv(diaumpire_path_list[expIdx], sep="\t", index_col=0, na_values=[0, "", "NA"], header=0, usecols=msfraggerdia_cols[expIdx])
    diann = pd.read_csv(diann_path_list[expIdx], sep="\t", index_col=0, na_values=[0, "", "NA"], header=0, usecols=msfraggerdia_cols[expIdx])
    maxdia = pd.read_csv(maxdia_path_list[expIdx], sep="\t", index_col=0, na_values=[0, "", "NaN"], header=0, usecols=["Protein IDs", "Intensity 1", "Intensity 2", "Intensity 3", "Reverse", "Potential contaminant"])
    msfraggerdia = pd.read_csv(msfraggerdia_path_list[expIdx], sep="\t", index_col=0, na_values=[0, "", "NA"], header=0, usecols=msfraggerdia_cols[expIdx])

    # remove rows with missing values to calculate CVs
    diaumpire.dropna(how="any", inplace=True)
    diann.dropna(how="any", inplace=True)

    maxdia = maxdia.loc[pd.isna(maxdia["Reverse"]) & pd.isna(maxdia["Potential contaminant"])]
    maxdia.drop(labels=["Reverse", "Potential contaminant"], axis=1, inplace=True)
    maxdia.dropna(how="any", inplace=True)

    msfraggerdia.dropna(how="any", inplace=True)

    diaumpire = diaumpire.groupby(level=0).max()
    diann = diann.groupby(level=0).max()
    maxdia = maxdia.groupby(level=0).max()
    msfraggerdia = msfraggerdia.groupby(level=0).max()

    total_counts = pd.DataFrame({"FP-DIAU": len(diaumpire.index), "DIA-NN\nlib-free": len(diann.index), "MaxDIA": len(maxdia.index), "FP-MSF": len(msfraggerdia.index)}, index=[0])

    diaumpire_cv = np.nanstd(diaumpire, 1) * 100 / np.nanmean(diaumpire, 1)
    diann_cv = np.nanstd(diann, 1) * 100 / np.nanmean(diann, 1)
    maxdia_cv = np.nanstd(maxdia, 1) * 100 / np.nanmean(maxdia, 1)
    msfraggerdia_cv = np.nanstd(msfraggerdia, 1) * 100 / np.nanmean(msfraggerdia, 1)

    total_counts2 = pd.DataFrame({"FP-DIAU": sum(diaumpire_cv < 20), "DIA-NN\nlib-free": sum(diann_cv < 20), "MaxDIA": len(maxdia_cv < 20), "FP-MSF": sum(msfraggerdia_cv < 20)}, index=[0])

    sns_plot = sns.barplot(data=total_counts, palette=[palette_dict["diaumpire"], palette_dict["diann"], palette_dict["maxdia"], palette_dict["msfraggerdia"]], alpha=0.6)
    sns_plot = sns.barplot(data=total_counts2, palette=[palette_dict["diaumpire"], palette_dict["diann"], palette_dict["maxdia"], palette_dict["msfraggerdia"]])
    sns_plot.set(xlabel=None)
    sns_plot.set_ylabel("quantified proteins", fontsize=13)
    sns_plot.set_xticklabels(sns_plot.get_xticklabels(), fontsize=11)
    sns_plot.set_title(experiment_list[expIdx], fontsize=15)
    sns_plot.figure.savefig("low_input_cell_protein_" + experiment_list[expIdx] + ".pdf")
    plt.figure()
