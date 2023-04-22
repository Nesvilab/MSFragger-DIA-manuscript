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

maxdia_path_list = [r"single-cell_1cell_10cells12xLib\maxdia\combined\txt\proteinGroups.txt",
                    r"single-cell_100cells_1ug8xLib\maxdia\combined\txt\proteinGroups.txt"]

experiment_list = ["1 cell", "117 cells"]
msfraggerdia_cols = [["Unnamed: 0", "20201112_MEC1-MF-1cell-Eclipse-DIA-01", "20201112_MEC1-MF-1cell-Eclipse-DIA-02", "20201112_MEC1-MF-1cell-Eclipse-DIA-03"],
                     ["Unnamed: 0", "20201112_MEC1-MF-100cells-Eclipse-DIA-01", "20201112_MEC1-MF-100cells-Eclipse-DIA-02", "20201112_MEC1-MF-100cells-Eclipse-DIA-03"]]

for expIdx in range(len(experiment_list)):
    maxdia_1 = pd.read_csv(maxdia_path_list[expIdx], sep="\t", index_col=0, na_values=[0, "", "NaN"], header=0, usecols=["Protein IDs", "Intensity 1", "Intensity 2", "Intensity 3", "Reverse", "Potential contaminant"])
    maxdia_2 = pd.read_csv(maxdia_path_list[expIdx], sep="\t", index_col=0, na_values=[0, "", "NaN"], header=0, usecols=["Protein IDs", "LFQ intensity 1", "LFQ intensity 2", "LFQ intensity 3", "Reverse", "Potential contaminant"])

    maxdia_1 = maxdia_1.loc[pd.isna(maxdia_1["Reverse"]) & pd.isna(maxdia_1["Potential contaminant"])]
    maxdia_1.drop(labels=["Reverse", "Potential contaminant"], axis=1, inplace=True)
    maxdia_1.dropna(how="all", inplace=True)

    maxdia_2 = maxdia_2.loc[pd.isna(maxdia_2["Reverse"]) & pd.isna(maxdia_2["Potential contaminant"])]
    maxdia_2.drop(labels=["Reverse", "Potential contaminant"], axis=1, inplace=True)
    maxdia_2.dropna(how="all", inplace=True)

    maxdia_1 = maxdia_1.groupby(level=0).max()
    maxdia_2 = maxdia_2.groupby(level=0).max()

    maxdia_count_1 = pd.DataFrame({"count": maxdia_1.count(axis=0), "tool": "MaxDIA intensity"})
    maxdia_count_2 = pd.DataFrame({"count": maxdia_2.count(axis=0), "tool": "MaxDIA MaxLFQ intensity"})

    total_counts = pd.DataFrame({"MaxDIA intensity": len(maxdia_1.index), "MaxDIA MaxLFQ intensity": len(maxdia_2.index)}, index=[0])

    run_counts = pd.concat([maxdia_count_1, maxdia_count_2])

    sns_plot = sns.barplot(data=total_counts, palette=[palette_dict["maxdia"], palette_dict["maxdia"]])
    sns_plot = sns.swarmplot(x="tool", y="count", data=run_counts, color="white")
    sns_plot.set(xlabel=None)
    sns_plot.set_ylabel("quantified proteins", fontsize=13)
    sns_plot.set_xticklabels(sns_plot.get_xticklabels(), fontsize=10)
    sns_plot.set_title(experiment_list[expIdx], fontsize=15)
    sns_plot.figure.savefig("single_cell_protein_" + experiment_list[expIdx] + "_S.pdf")
    plt.figure()
