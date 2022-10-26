import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


sns.set(font="Arial")
sns.set_theme(style="ticks")


# ccrcc 20 desktop
diaumpire = [254, 17, 23, 10]  # dia-umpire, msfragger, msbooster+percolator+filtering+easypqp, dia-nn
diann = [30, 332-30]  # prediction, others
maxdia = [1255]
msfraggerdia = [65, 50, 10]  # msfragger, msbooster+percolator+filtering+easypqp, dia-nn

xx = pd.DataFrame({"tool":  ["FP-DIAU", "DIA-NN", "MaxDIA", "FP-MSF"],
                   "spectral library prediction": [0, diann[0], 0, 0],
                   "identification and quantification": [0, diann[1], 0, 0],
                   "pseudo-MS/MS generation": [diaumpire[0], 0, 0, 0],
                   "database searching": [diaumpire[1], 0, 0, msfraggerdia[0]],
                   "rescoring and FDR filtering": [diaumpire[2], 0, 0, msfraggerdia[1]],
                   "quantification": [diaumpire[3], 0, 0, msfraggerdia[2]],
                   "MaxDIA": [0, 0, maxdia[0], 0]
                   })

fig = xx.set_index("tool").plot(kind="bar", stacked=True, xlabel="", ylabel="run time (minutes)", linewidth=0).get_figure()
plt.xticks(rotation=0)
fig.savefig("runtime_ccrcc_20_desktop.pdf")
plt.figure()


# ccrcc 20 dev02
diaumpire = [292, 15, 34, 17]  # dia-umpire, msfragger, msbooster+percolator+filtering+easypqp, dia-nn
diann = [10, 194]  # prediction, others
maxdia = [550]
msfraggerdia = [43, 47, 18]  # msfragger, msbooster+percolator+filtering+easypqp, dia-nn

xx = pd.DataFrame({"tool":  ["FP-DIAU", "DIA-NN", "MaxDIA", "FP-MSF"],
                   "spectral library prediction": [0, diann[0], 0, 0],
                   "identification and quantification": [0, diann[1], 0, 0],
                   "pseudo-MS/MS generation": [diaumpire[0], 0, 0, 0],
                   "database searching": [diaumpire[1], 0, 0, msfraggerdia[0]],
                   "rescoring and FDR filtering": [diaumpire[2], 0, 0, msfraggerdia[1]],
                   "quantification": [diaumpire[3], 0, 0, msfraggerdia[2]],
                   "MaxDIA": [0, 0, maxdia[0], 0]
                   })

fig = xx.set_index("tool").plot(kind="bar", stacked=True, xlabel="", ylabel="run time (minutes)", linewidth=0).get_figure()
plt.xticks(rotation=0)
fig.savefig("runtime_ccrcc_20_dev02.pdf")
plt.figure()


# melanoma desktop
diaumpire = [132, 13, 13, 5]  # dia-umpire, msfragger, msbooster+percolator+filtering+easypqp, dia-nn
diann = [1213, 2558-1213]  # prediction, others
msfraggerdia = [141, 33, 5]  # msfragger, msbooster+percolator+filtering+easypqp, dia-nn

xx = pd.DataFrame({"tool":  ["FP-DIAU", "DIA-NN", "FP-MSF"],
                   "spectral library prediction": [0, diann[0], 0],
                   "identification and quantification": [0, diann[1], 0],
                   "pseudo-MS/MS generation": [diaumpire[0], 0, 0],
                   "database searching using MSFragger": [diaumpire[1], 0, msfraggerdia[0]],
                   "rescoring and FDR filtering": [diaumpire[2], 0, msfraggerdia[1]],
                   "quantification": [diaumpire[3], 0, msfraggerdia[2]],
                   })

fig = xx.set_index("tool").plot(kind="bar", stacked=True, xlabel="", ylabel="run time (minutes)", linewidth=0).get_figure()
plt.xticks(rotation=0)
fig.savefig("runtime_melanoma_phospho_desktop.pdf")


# melanoma dev02
diaumpire = [157, 11, 14, 11]  # dia-umpire, msfragger, msbooster+percolator+filtering+easypqp, dia-nn
diann = [109, 1143-109]  # prediction, others
msfraggerdia = [144, 36, 13]  # msfragger, msbooster+percolator+filtering+easypqp, dia-nn

xx = pd.DataFrame({"tool":  ["FP-DIAU", "DIA-NN", "FP-MSF"],
                   "spectral library prediction": [0, diann[0], 0],
                   "identification and quantification": [0, diann[1], 0],
                   "pseudo-MS/MS generation": [diaumpire[0], 0, 0],
                   "database searching using MSFragger": [diaumpire[1], 0, msfraggerdia[0]],
                   "rescoring and FDR filtering": [diaumpire[2], 0, msfraggerdia[1]],
                   "quantification": [diaumpire[3], 0, msfraggerdia[2]],
                   })

fig = xx.set_index("tool").plot(kind="bar", stacked=True, xlabel="", ylabel="run time (minutes)", linewidth=0).get_figure()
plt.xticks(rotation=0)
fig.savefig("runtime_melanoma_phospho_dev02.pdf")

