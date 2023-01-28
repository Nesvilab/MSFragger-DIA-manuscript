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

diaumpire_path = r"G:\Dropbox\papers_Fengchao\msfragger_dia\script\results\ccRCC\diaumpire\gene_maxlfq.tsv"
diaumpiredda_path = r"G:\Dropbox\papers_Fengchao\msfragger_dia\script\results\ccRCC\diaumpiredda\gene_maxlfq.tsv"
diann_path = r"G:\Dropbox\papers_Fengchao\msfragger_dia\script\results\ccRCC\diann\gene_maxlfq.tsv"
maxdia_path = r"G:\Dropbox\papers_Fengchao\msfragger_dia\script\results\ccRCC\maxdia\combined\txt\proteinGroups.txt"
msfraggerdia_path = r"G:\Dropbox\papers_Fengchao\msfragger_dia\script\results\ccRCC\msfraggerdia\gene_maxlfq.tsv"
msfraggerdiadda_path = r"G:\Dropbox\papers_Fengchao\msfragger_dia\script\results\ccRCC\msfraggerdiadda\gene_maxlfq.tsv"

diaumpire = pd.read_csv(diaumpire_path, sep="\t", index_col=0, na_values=[0, "", "NA"], header=0)
diaumpiredda = pd.read_csv(diaumpiredda_path, sep="\t", index_col=0, na_values=[0, "", "NA"], header=0)
diann = pd.read_csv(diann_path, sep="\t", index_col=0, na_values=[0, "", "NA"], header=0)
maxdia = pd.read_csv(maxdia_path, sep="\t", index_col=0, na_values=[0, "", "NaN"], header=0, usecols=["Protein IDs", "Intensity L-00004_", "Intensity L-00004_NA", "Intensity L-00010_", "Intensity L-00010_NA", "Intensity L-00011_", "Intensity L-00011_NA", "Intensity L-00026_", "Intensity L-00026_NA", "Intensity L-00079_", "Intensity L-00079_NA", "Intensity L-00088_", "Intensity L-00088_NA", "Intensity L-00096_", "Intensity L-00096_NA", "Intensity L-00097_", "Intensity L-00097_NA", "Intensity L-00103_", "Intensity L-00103_NA", "Intensity L-00183_", "Intensity L-00183_NA", "Intensity L-00359_", "Intensity L-00360_", "Intensity L-00369_", "Intensity L-00369_NA", "Intensity L-00416_", "Intensity L-00416_NA", "Intensity L-00418_", "Intensity L-00418_NA", "Intensity L-00447_", "Intensity L-00447_NA", "Intensity L-00448_", "Intensity L-00448_NA", "Intensity L-00561_", "Intensity L-00561_NA", "Intensity L-00581_", "Intensity L-00581_NA", "Intensity L-00583_", "Intensity L-00583_NA", "Intensity L-00606_", "Intensity L-00606_NA", "Intensity L-00607_", "Intensity L-00607_NA", "Intensity L-00610_", "Intensity L-00765_", "Intensity L-00766_", "Intensity L-00790_", "Intensity L-00791_", "Intensity L-00791_NA", "Intensity L-00792_", "Intensity L-00796_", "Intensity L-00799_", "Intensity L-00800_", "Intensity L-00812_", "Intensity L-00813_", "Intensity L-00814_", "Intensity L-00814_NA", "Intensity L-00817_", "Intensity L-00902_", "Intensity L-00902_NA", "Intensity L-00907_", "Intensity L-00907_NA", "Intensity L-00908_", "Intensity L-00908_NA", "Intensity L-00910_", "Intensity L-01286_", "Intensity L-01286_NA", "Intensity L-01287_", "Intensity L-01287_NA", "Intensity L-01288_", "Intensity L-01302_", "Intensity L-01302_NA", "Intensity L-01313_", "Intensity L-01313_NA", "Intensity L-01352_", "Intensity L-01553_", "Intensity L-01557_", "Intensity L-01560_", "Intensity L-01603_", "Intensity L-01603_NA", "Intensity L-01607_", "Intensity L-01607_NA", "Intensity L-01836_", "Intensity L-01836_NA", "Intensity L-01861_", "Intensity L-01861_NA", "Intensity L-01882_", "Intensity L-01882_NA", "Intensity L-01885_", "Intensity L-01885_NA", "Intensity N-00148_", "Intensity N-00148_NA", "Intensity N-00149_", "Intensity N-00149_NA", "Intensity N-00150_", "Intensity N-00150_NA", "Intensity N-00154_", "Intensity N-00168_", "Intensity N-00168_NA", "Intensity N-00177_", "Intensity N-00177_NA", "Intensity N-00194_", "Intensity N-00194_NA", "Intensity N-00242_", "Intensity N-00242_NA", "Intensity N-00244_", "Intensity N-00244_NA", "Intensity N-00246_", "Intensity N-00246_NA", "Intensity N-00305_", "Intensity N-00310_", "Intensity N-00310_NA", "Intensity N-00312_", "Intensity N-00312_NA", "Intensity N-00313_", "Intensity N-00314_", "Intensity N-00314_NA", "Intensity N-00315_", "Intensity N-00317_", "Intensity N-00317_NA", "Intensity N-00320_", "Intensity N-00320_NA", "Intensity N-00380_", "Intensity N-00390_", "Intensity N-00390_NA", "Intensity N-00435_", "Intensity N-00435_NA", "Intensity N-00437_", "Intensity N-00437_NA", "Intensity N-00491_", "Intensity N-00491_NA", "Intensity N-00492_", "Intensity N-00492_NA", "Intensity N-00494_", "Intensity N-00494_NA", "Intensity N-00495_", "Intensity N-00495_NA", "Intensity N-00573_", "Intensity N-00573_NA", "Intensity N-00577_", "Intensity N-00577_NA", "Intensity N-00646_", "Intensity N-00646_NA", "Intensity N-00733_", "Intensity N-00733_NA", "Intensity N-00831_", "Intensity N-00831_NA", "Intensity N-00832_", "Intensity N-00834_", "Intensity N-00834_NA", "Intensity N-00852_", "Intensity N-00852_NA", "Intensity N-00953_", "Intensity N-00953_NA", "Intensity N-01175_", "Intensity N-01175_NA", "Intensity N-01176_", "Intensity N-01176_NA", "Intensity N-01178_", "Intensity N-01178_NA", "Intensity N-01179_", "Intensity N-01179_NA", "Intensity N-01180_", "Intensity N-01200_", "Intensity N-01200_NA", "Intensity N-01213_", "Intensity N-01214_", "Intensity N-01214_NA", "Intensity N-01220_", "Intensity N-01220_NA", "Intensity N-01261_", "Intensity N-01261_NA", "Intensity N-01361_", "Intensity N-01361_NA", "Intensity N-01522_", "Intensity N-01522_NA", "Intensity N-01524_", "Intensity N-01524_NA", "Intensity N-01646_", "Intensity N-01646_NA", "Intensity N-01648_", "Intensity N-01648_NA", "Intensity N-01649_", "Intensity N-01649_NA", "Intensity N-01651_", "Intensity N-01651_NA", "Intensity N-01808_", "Intensity N-01808_NA", "Reverse", "Potential contaminant"], low_memory=False)
msfraggerdia = pd.read_csv(msfraggerdia_path, sep="\t", index_col=0, na_values=[0, "", "NA"], header=0)
msfraggerdiadda = pd.read_csv(msfraggerdiadda_path, sep="\t", index_col=0, na_values=[0, "", "NA"], header=0)

diaumpire.dropna(thresh=175*0.5, inplace=True)
diaumpiredda.dropna(thresh=175*0.5, inplace=True)
diann.dropna(thresh=175*0.5, inplace=True)

maxdia = maxdia.loc[pd.isna(maxdia["Reverse"]) & pd.isna(maxdia["Potential contaminant"])]
maxdia.drop(labels=["Reverse", "Potential contaminant"], axis=1, inplace=True)
maxdia.dropna(thresh=175*0.5, inplace=True)

msfraggerdia.dropna(thresh=175*0.5, inplace=True)
msfraggerdiadda.dropna(thresh=175*0.5, inplace=True)

diaumpire = diaumpire.groupby(level=0).max()
diaumpiredda = diaumpiredda.groupby(level=0).max()
diann = diann.groupby(level=0).max()
maxdia = maxdia.groupby(level=0).max()
msfraggerdia = msfraggerdia.groupby(level=0).max()
msfraggerdiadda = msfraggerdiadda.groupby(level=0).max()

diaumpire_count = pd.DataFrame({"count": diaumpire.count(axis=0), "tool": "FP-DIAU"})
diaumpiredda_count = pd.DataFrame({"count": diaumpiredda.count(axis=0), "tool": "FP-DIAU\nhybrid"})
diann_count = pd.DataFrame({"count": diann.count(axis=0), "tool": "DIA-NN"})
maxdia_count = pd.DataFrame({"count": maxdia.count(axis=0), "tool": "MaxDIA"})
msfraggerdia_count = pd.DataFrame({"count": msfraggerdia.count(axis=0), "tool": "FP-MSF"})
msfraggerdiadda_count = pd.DataFrame({"count": msfraggerdiadda.count(axis=0), "tool": "FP-MSF\nhybrid"})

total_counts = pd.DataFrame({"FP-DIAU": len(diaumpire.index), "FP-DIAU\nhybrid": len(diaumpiredda.index), "DIA-NN": len(diann.index), "MaxDIA": len(maxdia.index), "FP-MSF": len(msfraggerdia.index), "FP-MSF\nhybrid": len(msfraggerdiadda.index)}, index=[0])

run_counts = pd.concat([diaumpire_count, diaumpiredda_count, diann_count, maxdia_count, msfraggerdia_count, msfraggerdiadda_count])

sns_plot = sns.barplot(data=total_counts, palette=[palette_dict["diaumpire"], palette_dict["diaumpire hybrid"], palette_dict["diann"], palette_dict["maxdia"], palette_dict["msfraggerdia"], palette_dict["msfragger hybrid"]])
sns_plot = sns.swarmplot(x="tool", y="count", data=run_counts, color="white", size=1.3)
sns_plot.set(xlabel=None, ylabel="quantified genes")
sns_plot.set_xticklabels(sns_plot.get_xticklabels(), fontsize=11)
sns_plot.figure.savefig("ccRCC.pdf")
