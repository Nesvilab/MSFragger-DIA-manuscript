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
    "diann": sns.color_palette()[3],
    "maxdia": sns.color_palette()[2],
    "msfraggerdia": sns.color_palette()[0],
    "msfragger hybrid": sns.color_palette()[1]
}

# plot a
diann_path = r"G:\Dropbox\papers_Fengchao\msfragger_dia\script\results\2020-Yeast\diann\modified_sequence_maxlfq.tsv"
maxdia_path = r"G:\Dropbox\papers_Fengchao\msfragger_dia\script\results\2020-Yeast\maxdia\combined\txt\modificationSpecificPeptides.txt"
msfraggerdia_path = r"G:\Dropbox\papers_Fengchao\msfragger_dia\script\results\2020-Yeast\msfraggerdia\modified_sequence_maxlfq.tsv"

diann = pd.read_csv(diann_path, sep="\t", index_col=0, na_values=[0, "", "NA"], header=0, usecols=["Unnamed: 0", "20190206_LUM1_CPBA_EASY04_060_30_SA_90mingrad_80B_DIA_400_1000_8mzol_15k_20IIT_4e5agc_1633-01_01", "20190206_LUM1_CPBA_EASY04_060_30_SA_90mingrad_80B_DIA_400_1000_8mzol_15k_20IIT_4e5agc_1633-01_02", "20190206_LUM1_CPBA_EASY04_060_30_SA_90mingrad_80B_DIA_400_1000_8mzol_15k_20IIT_4e5agc_1633-01_04", "20190206_LUM1_CPBA_EASY04_060_30_SA_90mingrad_80B_DIA_400_1000_8mzol_15k_20IIT_4e5agc_1633-01_05"])

maxdia = pd.read_csv(maxdia_path, sep="\t", index_col=0, na_values=[0, "", "NaN"], header=0, usecols=["Sequence", "Intensity 20190206_LUM1_CPBA_EASY04_060_30_SA_90mingrad_80B_DIA_400_1000_8mzol_15k_20IIT_4e5agc_1633-01_01", "Intensity 20190206_LUM1_CPBA_EASY04_060_30_SA_90mingrad_80B_DIA_400_1000_8mzol_15k_20IIT_4e5agc_1633-01_02", "Intensity 20190206_LUM1_CPBA_EASY04_060_30_SA_90mingrad_80B_DIA_400_1000_8mzol_15k_20IIT_4e5agc_1633-01_04", "Intensity 20190206_LUM1_CPBA_EASY04_060_30_SA_90mingrad_80B_DIA_400_1000_8mzol_15k_20IIT_4e5agc_1633-01_05", "Reverse", "Potential contaminant"])

msfraggerdia = pd.read_csv(msfraggerdia_path, sep="\t", index_col=0, na_values=[0, "", "NA"], header=0, usecols=["Unnamed: 0", "20190206_LUM1_CPBA_EASY04_060_30_SA_90mingrad_80B_DIA_400_1000_8mzol_15k_20IIT_4e5agc_1633-01_01", "20190206_LUM1_CPBA_EASY04_060_30_SA_90mingrad_80B_DIA_400_1000_8mzol_15k_20IIT_4e5agc_1633-01_02", "20190206_LUM1_CPBA_EASY04_060_30_SA_90mingrad_80B_DIA_400_1000_8mzol_15k_20IIT_4e5agc_1633-01_04", "20190206_LUM1_CPBA_EASY04_060_30_SA_90mingrad_80B_DIA_400_1000_8mzol_15k_20IIT_4e5agc_1633-01_05"])

diann.dropna(how="all", inplace=True)

maxdia = maxdia.loc[pd.isna(maxdia["Reverse"]) & pd.isna(maxdia["Potential contaminant"])]
maxdia.drop(labels=["Reverse", "Potential contaminant"], axis=1, inplace=True)
maxdia.dropna(how="all", inplace=True)

msfraggerdia.dropna(how="all", inplace=True)

publication_count = pd.DataFrame({"count": [44676, 45403, 41987, 41511], "tool": "Searle et al. 2020"})
diann_count = pd.DataFrame({"count": diann.count(axis=0), "tool": "DIA-NN"})
maxdia_count = pd.DataFrame({"count": maxdia.count(axis=0), "tool": "MaxDIA"})
msfraggerdia_count = pd.DataFrame({"count": msfraggerdia.count(axis=0), "tool": "FP-MSF"})

average_counts = pd.DataFrame({"Searle et al. 2020": publication_count["count"].mean(), "DIA-NN": diann_count["count"].mean(), "FP-MSF": msfraggerdia_count["count"].mean()}, index=[0])

run_counts = pd.concat([publication_count, diann_count, msfraggerdia_count])

sns_plot = sns.barplot(data=average_counts, palette=[palette_dict["publication"], palette_dict["diann"], palette_dict["msfraggerdia"]])
sns_plot = sns.swarmplot(x="tool", y="count", data=run_counts, color="white", edgecolor="black", size=10, linewidth=2)
sns_plot.set(xlabel=None, ylabel="quantified peptides")
sns_plot.set_xticklabels(sns_plot.get_xticklabels())
sns_plot.figure.savefig("2020-Yeast_a.pdf")
plt.figure()


# plot a S
average_counts = pd.DataFrame({"Searle et al. 2020": publication_count["count"].mean(), "DIA-NN": diann_count["count"].mean(), "MaxDIA": maxdia_count["count"].mean(), "FP-MSF": msfraggerdia_count["count"].mean()}, index=[0])

run_counts = pd.concat([publication_count, diann_count, maxdia_count, msfraggerdia_count])

sns_plot = sns.barplot(data=average_counts, palette=[palette_dict["publication"], palette_dict["diann"], palette_dict["maxdia"], palette_dict["msfraggerdia"]])
sns_plot = sns.swarmplot(x="tool", y="count", data=run_counts, color="white", edgecolor="black", size=10, linewidth=2)
sns_plot.set(xlabel=None, ylabel="quantified peptides")
sns_plot.set_xticklabels(sns_plot.get_xticklabels())
sns_plot.figure.savefig("2020-Yeast_a_S.pdf")
plt.figure()


# plot b
xx = diann.dropna(axis=0, thresh=2)
diann_cv = pd.DataFrame({"cv": np.nanstd(xx, 1) * 100 / np.nanmean(xx, 1), "tool": "DIA-NN"})
xx = maxdia.dropna(axis=0, thresh=2)
maxdia_cv = pd.DataFrame({"cv": np.nanstd(xx, 1) * 100 / np.nanmean(xx, 1), "tool": "MaxDIA"})
xx = msfraggerdia.dropna(axis=0, thresh=2)
msfraggerdia_cv = pd.DataFrame({"cv": np.nanstd(xx, 1) * 100 / np.nanmean(xx, 1), "tool": "FP-MSF"})

x = pd.concat([diann_cv, msfraggerdia_cv])

sns_plot = sns.violinplot(x="tool", y="cv", data=x, palette=[palette_dict["diann"], palette_dict["msfraggerdia"]])
sns_plot.set(xlabel=None, ylabel="coefficient of variation (%)")
sns_plot.figure.savefig("2020-Yeast_b.pdf")
plt.figure()


# plot b S
x = pd.concat([diann_cv, maxdia_cv, msfraggerdia_cv])

sns_plot = sns.violinplot(x="tool", y="cv", data=x, palette=[palette_dict["diann"], palette_dict["maxdia"], palette_dict["msfraggerdia"]])
sns_plot.set(xlabel=None, ylabel="coefficient of variation (%)")
sns_plot.figure.savefig("2020-Yeast_b_S.pdf")
