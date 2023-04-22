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

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from upsetplot import from_contents, UpSet

sns.set(font="Arial")
sns.set_theme(style="ticks")


def translate_spectronaut_peptide(aa):
    aa.index = aa.index.str.\
        replace("[Acetyl (Protein N-term)]", "(UniMod:1)", regex=False).str.\
        replace("[Oxidation (M)]", "(UniMod:35)", regex=False).str.\
        replace("[Carbamidomethyl (C)]", "(UniMod:4)", regex=False).str.\
        replace("[Phospho (STY)]", "(UniMod:21)", regex=False).str.\
        replace("_", "", regex=False).str. \
        replace(".", "", regex=False)
    return aa


def translate_spectronaut_protein(aa):
    aa.index = aa.index.str.replace(";.+", "", regex=True)
    return aa


os.chdir(r"G:\Dropbox\papers_Fengchao\msfragger_dia\script\results\benchmark")

# precursor level
spectronaut_14_path = r"spectronaut\14\precursor_maxlfq.tsv"
spectronaut_17_path = r"spectronaut\17\precursor_maxlfq.tsv"
diann_path = r"diann\precursor_maxlfq.tsv"
msfraggerdia_path = r"msfraggerdia\precursor_maxlfq.tsv"
msfraggerdia_hybrid_path = r"msfraggerdiadda\precursor_maxlfq.tsv"

spectronaut_14 = pd.read_csv(spectronaut_14_path, sep="\t", index_col=0, na_values="NA", header=0)
spectronaut_17 = pd.read_csv(spectronaut_17_path, sep="\t", index_col=0, na_values="NA", header=0)
diann = pd.read_csv(diann_path, sep="\t", index_col=0, na_values="NA", header=0)
msfraggerdia = pd.read_csv(msfraggerdia_path, sep="\t", index_col=0, na_values="NA", header=0)
msfraggerdia_hybrid = pd.read_csv(msfraggerdia_hybrid_path, sep="\t", index_col=0, na_values="NA", header=0)

spectronaut_14.dropna(how="all", inplace=True)
spectronaut_17.dropna(how="all", inplace=True)
diann.dropna(how="all", inplace=True)
msfraggerdia.dropna(how="all", inplace=True)
msfraggerdia_hybrid.dropna(how="all", inplace=True)

spectronaut_14 = translate_spectronaut_peptide(spectronaut_14)
spectronaut_17 = translate_spectronaut_peptide(spectronaut_17)

spectronaut_14.drop_duplicates(keep="first", inplace=True)
spectronaut_17.drop_duplicates(keep="first", inplace=True)
diann.drop_duplicates(keep="first", inplace=True)
msfraggerdia.drop_duplicates(keep="first", inplace=True)
msfraggerdia_hybrid.drop_duplicates(keep="first", inplace=True)

tt = from_contents({"Spectronaut 14": spectronaut_14.index.tolist(), "Spectronaut 17": spectronaut_17.index.tolist(), "DIA-NN lib-free": diann.index.tolist(), "FP-MSF": msfraggerdia.index.tolist(), "FP-MSF hybrid": msfraggerdia_hybrid.index.tolist()})

UpSet(tt, subset_size="count").plot(y_label="precursors")
plt.savefig("benchmark_upset_precursor.pdf")
