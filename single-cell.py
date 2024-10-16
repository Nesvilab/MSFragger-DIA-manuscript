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

diaumpiredda_path_list = [r"single-cell_1cell_10cells12xLib\diaumpiredda\protein_maxlfq.tsv", r"single-cell_100cells_1ug8xLib\diaumpiredda\protein_maxlfq.tsv"]
diann_path_list = [r"single-cell_1cell_10cells12xLib\diann\protein_maxlfq.tsv", r"single-cell_100cells_1ug8xLib\diann\protein_maxlfq.tsv"]
maxdia_path_list = [r"single-cell_1cell_10cells12xLib\maxdia\combined\txt\proteinGroups.txt", r"single-cell_100cells_1ug8xLib\maxdia\combined\txt\proteinGroups.txt"]
msfraggerdia_path_list = [r"single-cell_1cell_10cells12xLib\msfraggerdia\protein_maxlfq.tsv", r"single-cell_100cells_1ug8xLib\msfraggerdia\protein_maxlfq.tsv"]
msfraggerdiadda_path_list = [r"single-cell_1cell_10cells12xLib\msfraggerdiadda\protein_maxlfq.tsv", r"single-cell_100cells_1ug8xLib\msfraggerdiadda\protein_maxlfq.tsv"]

experiment_list = ["1 cell", "117 cells"]
msfraggerdia_cols = [["Unnamed: 0", "20201112_MEC1-MF-1cell-Eclipse-DIA-01", "20201112_MEC1-MF-1cell-Eclipse-DIA-02", "20201112_MEC1-MF-1cell-Eclipse-DIA-03"],
                     ["Unnamed: 0", "20201112_MEC1-MF-100cells-Eclipse-DIA-01", "20201112_MEC1-MF-100cells-Eclipse-DIA-02", "20201112_MEC1-MF-100cells-Eclipse-DIA-03"]]

for expIdx in range(len(experiment_list)):
    diaumpiredda = pd.read_csv(diaumpiredda_path_list[expIdx], sep="\t", index_col=0, na_values="NA", header=0, usecols=msfraggerdia_cols[expIdx])
    diann = pd.read_csv(diann_path_list[expIdx], sep="\t", index_col=0, na_values="NA", header=0, usecols=msfraggerdia_cols[expIdx])
    maxdia = pd.read_csv(maxdia_path_list[expIdx], sep="\t", index_col=0, na_values=[0, "", "NaN"], header=0, usecols=["Protein IDs", "Intensity 1", "Intensity 2", "Intensity 3", "Reverse", "Potential contaminant"])
    msfraggerdia = pd.read_csv(msfraggerdia_path_list[expIdx], sep="\t", index_col=0, na_values="NA", header=0, usecols=msfraggerdia_cols[expIdx])
    msfraggerdiadda = pd.read_csv(msfraggerdiadda_path_list[expIdx], sep="\t", index_col=0, na_values="NA", header=0, usecols=msfraggerdia_cols[expIdx])

    diaumpiredda.dropna(how="all", inplace=True)
    diann.dropna(how="all", inplace=True)

    maxdia = maxdia.loc[pd.isna(maxdia["Reverse"]) & pd.isna(maxdia["Potential contaminant"])]
    maxdia.drop(labels=["Reverse", "Potential contaminant"], axis=1, inplace=True)
    maxdia.dropna(how="all", inplace=True)

    msfraggerdia.dropna(how="all", inplace=True)
    msfraggerdiadda.dropna(how="all", inplace=True)

    diaumpiredda = diaumpiredda.groupby(level=0).max()
    diann = diann.groupby(level=0).max()
    maxdia = maxdia.groupby(level=0).max()
    msfraggerdia = msfraggerdia.groupby(level=0).max()
    msfraggerdiadda = msfraggerdiadda.groupby(level=0).max()

    diaumpiredda_count = pd.DataFrame({"count": diaumpiredda.count(axis=0), "tool": "FP-DIAU\nhybrid"})
    diann_count = pd.DataFrame({"count": diann.count(axis=0), "tool": "DIA-NN\nlib-free"})
    maxdia_count = pd.DataFrame({"count": maxdia.count(axis=0), "tool": "MaxDIA"})
    msfraggerdia_count = pd.DataFrame({"count": msfraggerdia.count(axis=0), "tool": "FP-MSF"})
    msfraggerdiadda_count = pd.DataFrame({"count": msfraggerdiadda.count(axis=0), "tool": "FP-MSF\nhybrid"})

    total_counts = pd.DataFrame({"FP-DIAU\nhybrid": len(diaumpiredda.index), "DIA-NN\nlib-free": len(diann.index), "MaxDIA": len(maxdia.index), "FP-MSF": len(msfraggerdia.index), "FP-MSF\nhybrid": len(msfraggerdiadda.index)}, index=[0])

    run_counts = pd.concat([diaumpiredda_count, diann_count, maxdia_count, msfraggerdia_count, msfraggerdiadda_count])

    sns_plot = sns.barplot(data=total_counts, palette=[palette_dict["diaumpire hybrid"], palette_dict["diann"], palette_dict["maxdia"], palette_dict["msfraggerdia"], palette_dict["msfragger hybrid"]])
    sns_plot = sns.swarmplot(x="tool", y="count", data=run_counts, color="white")
    sns_plot.set(xlabel=None)
    sns_plot.set_ylabel("quantified proteins", fontsize=13)
    sns_plot.set_xticklabels(sns_plot.get_xticklabels(), fontsize=10)
    sns_plot.set_title(experiment_list[expIdx], fontsize=15)
    sns_plot.figure.savefig("single_cell_protein_" + experiment_list[expIdx] + ".pdf")
    plt.figure()
