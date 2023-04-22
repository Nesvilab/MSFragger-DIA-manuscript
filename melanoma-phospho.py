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
from upsetplot import from_contents, UpSet
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


def collapse_peptide(df):
    peptides = list(df.index.str.replace("\\(UniMod:[0-9]+\\)", "", regex=True))
    df["peptide"] = peptides
    return df.groupby("peptide").max()


os.chdir(r"G:\Dropbox\papers_Fengchao\msfragger_dia\script\results\melanoma-phospho")

diaumpire_path = r"diaumpire\modified_sequence_maxlfq.tsv"
diann_path = r"diann\modified_sequence_maxlfq.tsv"
msfraggerdia_path = r"msfraggerdia\modified_sequence_maxlfq.tsv"

diaumpire = pd.read_csv(diaumpire_path, sep="\t", index_col=0, na_values="NA", header=0)
diann = pd.read_csv(diann_path, sep="\t", index_col=0, na_values="NA", header=0)
msfraggerdia = pd.read_csv(msfraggerdia_path, sep="\t", index_col=0, na_values="NA", header=0)

diaumpire.dropna(thresh=6, inplace=True)
diann.dropna(thresh=6, inplace=True)
msfraggerdia.dropna(thresh=6, inplace=True)

diaumpire = diaumpire.groupby(level=0).max()
diann = diann.groupby(level=0).max()
msfraggerdia = msfraggerdia.groupby(level=0).max()

diaumpire = diaumpire.loc[diaumpire.index.str.contains("UniMod:21")]
diann = diann.loc[diann.index.str.contains("UniMod:21")]
msfraggerdia = msfraggerdia.loc[msfraggerdia.index.str.contains("UniMod:21")]

diaumpire = collapse_peptide(diaumpire)
diann = collapse_peptide(diann)
msfraggerdia = collapse_peptide(msfraggerdia)

ttt = from_contents({"DIA-Umpire": diaumpire.index.tolist(), "DIA-NN lib-free": diann.index.tolist(), "FP-MSF": msfraggerdia.index.tolist()})

UpSet(ttt, subset_size="count").plot(y_label="phosphosequences")
plt.savefig("melanoma_upset_peptide.pdf")
