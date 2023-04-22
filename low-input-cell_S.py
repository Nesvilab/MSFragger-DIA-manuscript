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

os.chdir(r"G:\Dropbox\papers_Fengchao\msfragger_dia\script\results")

maxdia_path_list = [r"low-input-cell_0.75ng_1.5ng\maxdia\combined\txt\proteinGroups.txt",
                    r"low-input-cell_7.5ng_1ug\maxdia\combined\txt\proteinGroups.txt"]

experiment_list = ["0.75 ng", "7.5 ng"]
msfraggerdia_cols = [["Unnamed: 0", "20210430_PC9-750pg_60K-IT118-DIA-IW10_01", "20210430_PC9-750pg_60K-IT118-DIA-IW10_02", "20210430_PC9-750pg_60K-IT118-DIA-IW10_03"],
                     ["Unnamed: 0", "20210430_PC9-7500pg_60K-IT118-DIA-IW10_01", "20210430_PC9-7500pg_60K-IT118-DIA-IW10_02", "20210430_PC9-7500pg_60K-IT118-DIA-IW10_03"]]

for expIdx in range(len(experiment_list)):
    maxdia_1 = pd.read_csv(maxdia_path_list[expIdx], sep="\t", index_col=0, na_values=[0, "", "NaN"], header=0, usecols=["Protein IDs", "Intensity 1", "Intensity 2", "Intensity 3", "Reverse", "Potential contaminant"])
    maxdia_2 = pd.read_csv(maxdia_path_list[expIdx], sep="\t", index_col=0, na_values=[0, "", "NaN"], header=0, usecols=["Protein IDs", "LFQ intensity 1", "LFQ intensity 2", "LFQ intensity 3", "Reverse", "Potential contaminant"])

    maxdia_1 = maxdia_1.loc[pd.isna(maxdia_1["Reverse"]) & pd.isna(maxdia_1["Potential contaminant"])]
    maxdia_1.drop(labels=["Reverse", "Potential contaminant"], axis=1, inplace=True)
    maxdia_1.dropna(how="any", inplace=True)

    maxdia_2 = maxdia_2.loc[pd.isna(maxdia_2["Reverse"]) & pd.isna(maxdia_2["Potential contaminant"])]
    maxdia_2.drop(labels=["Reverse", "Potential contaminant"], axis=1, inplace=True)
    maxdia_2.dropna(how="any", inplace=True)

    maxdia_1 = maxdia_1.groupby(level=0).max()
    maxdia_2 = maxdia_2.groupby(level=0).max()

    total_counts = pd.DataFrame({"MaxDIA intensity": len(maxdia_1.index), "MaxDIA MaxLFQ intensity": len(maxdia_2.index)}, index=[0])

    maxdia_cv_1 = np.nanstd(maxdia_1, 1) * 100 / np.nanmean(maxdia_1, 1)
    maxdia_cv_2 = np.nanstd(maxdia_2, 1) * 100 / np.nanmean(maxdia_2, 1)

    total_counts2 = pd.DataFrame({"MaxDIA intensity": len(maxdia_cv_1 < 20), "MaxDIA MaxLFQ intensity": len(maxdia_cv_2 < 20)}, index=[0])

    sns_plot = sns.barplot(data=total_counts, palette=[palette_dict["maxdia"], palette_dict["maxdia"]], alpha=0.6)
    sns_plot = sns.barplot(data=total_counts2, palette=[palette_dict["maxdia"], palette_dict["maxdia"]])
    sns_plot.set(xlabel=None)
    sns_plot.set_ylabel("quantified proteins", fontsize=13)
    sns_plot.set_xticklabels(sns_plot.get_xticklabels(), fontsize=11)
    sns_plot.set_title(experiment_list[expIdx], fontsize=15)
    sns_plot.figure.savefig("low_input_cell_protein_" + experiment_list[expIdx] + "_S.pdf")
    plt.figure()
