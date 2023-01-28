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
import numpy as np
import matplotlib.pyplot as plt

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


def collapse_peptide(df):
    peptides = list(df.index.str.replace("\\(UniMod:[0-9]+\\)", "", regex=True))
    df["peptide"] = peptides
    return df.groupby("peptide").max()


def collapse_spectronaut_peptide(df):
    peptides = list(df.index.str.replace("\\[[^\\[\\]]+\\]", "", regex=True))
    df["peptide"] = peptides
    return df.groupby("peptide").max()


diaumpire_path = r"G:\Dropbox\papers_Fengchao\msfragger_dia\script\results\melanoma-phospho\diaumpire\modified_sequence_maxlfq.tsv"
spectronaut_path = r"G:\Dropbox\papers_Fengchao\msfragger_dia\script\results\melanoma-phospho\spectronaut\20220105_161815_SPnew_Melanoma_directDIA_phosLocalized_20201020_Peptide_Report.tsv"
diann_path = r"G:\Dropbox\papers_Fengchao\msfragger_dia\script\results\melanoma-phospho\diann\modified_sequence_maxlfq.tsv"
msfraggerdia_path = r"G:\Dropbox\papers_Fengchao\msfragger_dia\script\results\melanoma-phospho\msfraggerdia\modified_sequence_maxlfq.tsv"

diaumpire = pd.read_csv(diaumpire_path, sep="\t", index_col=0, na_values=[0, "", "NA"], header=0)
spectronaut = pd.read_csv(spectronaut_path, sep="\t", index_col=0, na_values=[0, "", "Filtered"], header=0)
diann = pd.read_csv(diann_path, sep="\t", index_col=0, na_values=[0, "", "NA"], header=0)
msfraggerdia = pd.read_csv(msfraggerdia_path, sep="\t", index_col=0, na_values=[0, "", "NA"], header=0)

diaumpire.dropna(thresh=6, inplace=True)
spectronaut.dropna(thresh=6, inplace=True)
diann.dropna(thresh=6, inplace=True)
msfraggerdia.dropna(thresh=6, inplace=True)

diaumpire = diaumpire.groupby(level=0).max()
spectronaut = spectronaut.groupby(level=0).max()
diann = diann.groupby(level=0).max()
msfraggerdia = msfraggerdia.groupby(level=0).max()

diaumpire = diaumpire.loc[diaumpire.index.str.contains("UniMod:21")]
spectronaut = spectronaut.loc[spectronaut.index.str.contains("Phospho")]
diann = diann.loc[diann.index.str.contains("UniMod:21")]
msfraggerdia = msfraggerdia.loc[msfraggerdia.index.str.contains("UniMod:21")]

diaumpire = collapse_peptide(diaumpire)
spectronaut = collapse_spectronaut_peptide(spectronaut)
diann = collapse_peptide(diann)
msfraggerdia = collapse_peptide(msfraggerdia)

diaumpire_count = pd.DataFrame({"count": diaumpire.count(axis=0), "tool": "FP-DIAU"})
spectronaut_count = pd.DataFrame({"count": spectronaut.count(axis=0), "tool": "Spectronaut"})
diann_count = pd.DataFrame({"count": diann.count(axis=0), "tool": "DIA-NN"})
msfraggerdia_count = pd.DataFrame({"count": msfraggerdia.count(axis=0), "tool": "FP-MSF"})

total_counts = pd.DataFrame({"FP-DIAU": len(diaumpire.index), "Spectronaut": len(spectronaut.index), "DIA-NN": len(diann.index), "FP-MSF": len(msfraggerdia.index)}, index=[0])

run_counts = pd.concat([diaumpire_count, spectronaut_count, diann_count, msfraggerdia_count])

sns_plot = sns.barplot(data=total_counts, palette=[palette_dict["diaumpire"], palette_dict["spectronaut"], palette_dict["diann"], palette_dict["msfraggerdia"]])
sns_plot = sns.swarmplot(x="tool", y="count", data=run_counts, color="white", size=6)
sns_plot.set(xlabel=None, ylabel="quantified phosphosequences")
sns_plot.set_xticklabels(sns_plot.get_xticklabels(), fontsize=13)
sns_plot.figure.savefig("melanoma-phospho.pdf")
plt.figure()
