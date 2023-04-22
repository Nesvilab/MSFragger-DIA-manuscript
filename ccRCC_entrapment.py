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
import os


sns.set(font="Arial")
sns.set_theme(style="ticks")


def get_true_false_count(protein_column_list, protein_map):
    true_count = 0
    false_count = 0
    protein_list = [ss.split(";") for ss in protein_column_list]
    for proteins in protein_list:
        protein_type = -1
        # As long as there ia a Human protein in the protein group, it is a true protein group.
        for protein in proteins:
            if protein_map[protein].endswith("_HUMAN"):
                protein_type = 1
                break
        if protein_type == 1:
            true_count += 1
        elif protein_type == -1:
            false_count += 1
    return true_count, false_count


os.chdir(r"G:\Dropbox\papers_Fengchao\msfragger_dia\script\results\ccRCC_entrapment")

diann_path = r"diann\protein_maxlfq.tsv"
msfraggerdia_path = r"msfraggerdia\protein_maxlfq.tsv"
msfraggerdiadda_path = r"msfraggerdiadda\protein_maxlfq.tsv"
fasta_path = r"G:\Dropbox\papers_Fengchao\msfragger_dia\script\results\ccRCC_entrapment\2023-04-17-reviewed-UP000005640-UP000002311-UP000006548-UP000000625.fas"

human_protein_counts = 0
other_protein_counts = 0
protein_map = {}
for line in open(fasta_path).readlines():
    if line.startswith(">"):
        if "OS=Homo sapiens" in line:
            human_protein_counts += 1
        else:
            other_protein_counts += 1

        x = line[1::].strip().split(" ")[0].split("|")
        if len(x) > 2:
            protein_map[x[1]] = x[2]
        else:  # non-uniprot proteins such as iRT
            protein_map[line[1::].strip().split(" ")[0]] = ""

diann = pd.read_csv(diann_path, sep="\t", index_col=0, na_values="NA", header=0)
msfraggerdia = pd.read_csv(msfraggerdia_path, sep="\t", index_col=0, na_values="NA", header=0)
msfraggerdiadda = pd.read_csv(msfraggerdiadda_path, sep="\t", index_col=0, na_values="NA", header=0)

diann = diann.groupby(level=0).max()
msfraggerdia = msfraggerdia.groupby(level=0).max()
msfraggerdiadda = msfraggerdiadda.groupby(level=0).max()

diann_true_count, diann_false_count = get_true_false_count(diann.index.tolist(), protein_map)
msfraggerdia_true_count, msfraggerdia_false_count = get_true_false_count(msfraggerdia.index.tolist(), protein_map)
msfraggerdiadda_true_count, msfraggerdiadda_false_count = get_true_false_count(msfraggerdiadda.index.tolist(), protein_map)

print("human protein counts: " + str(human_protein_counts))
print("other protein counts: " + str(other_protein_counts))

print("diann true count: " + str(diann_true_count))
print("diann false count: " + str(diann_false_count))
print("diann FDP: " + str((diann_false_count * human_protein_counts) / (diann_true_count * other_protein_counts)))
print("msfraggerdia true count: " + str(msfraggerdia_true_count))
print("msfraggerdia false count: " + str(msfraggerdia_false_count))
print("msfraggerdia FDP: " + str((msfraggerdia_false_count * human_protein_counts) / (msfraggerdia_true_count * other_protein_counts)))
print("msfraggerdiadda true count: " + str(msfraggerdiadda_true_count))
print("msfraggerdiadda false count: " + str(msfraggerdiadda_false_count))
print("msfraggerdiadda FDP: " + str((msfraggerdiadda_false_count * human_protein_counts) / (msfraggerdiadda_true_count * other_protein_counts)))
