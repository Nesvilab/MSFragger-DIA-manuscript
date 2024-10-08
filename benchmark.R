################################################################################
############ Analysis of 92 Lymphnode + E.coli spike-in dataset ################
################################################################################
# by Klemens Froehlich 2022.02.14
# new DIA-NN and Spectronaut Version + MSFragger/DIA-NN integration
# modified by Fengchao Yu

rm(list = ls())

#load libraries
library(ggplot2)
library(tidyverse)
library(data.table)

theme_set(theme_classic())

setwd("G:\\Dropbox\\papers_Fengchao\\msfragger_dia\\script\\results\\benchmark\\")

# load annotation data
annotation <- read_tsv("G:\\Dropbox\\papers_Fengchao\\msfragger_dia\\script\\code\\MSFragger-DIA-manuscript\\uniprot.tab", col_select = c("Entry", "Organism"))

human_ecoli_overlap_peptide <- read.delim("human_ecoli_overlap_peptides.txt", header = FALSE, sep = "\t")

#### list library stats ####
lib1 <- read_tsv("msfraggerdia\\library.tsv", col_select = c(ProteinId, PeptideSequence, ModifiedPeptideSequence, PrecursorCharge))
lib2 <- lib1 %>%
  unique() %>%
  filter(!(PeptideSequence %in% human_ecoli_overlap_peptide$V1)) %>%
  mutate(PrecursorId = paste0(ModifiedPeptideSequence, PrecursorCharge)) %>%
  separate(col = ProteinId, into = "Entry", sep = ";", remove = TRUE, extra = "drop") %>%
  inner_join(annotation)
msfraggerdia_lib_precursors <- lib2 %>%
  group_by(Organism) %>%
  summarise(Precursors = n_distinct(PrecursorId))
msfraggerdia_lib_proteins <- lib2 %>%
  group_by(Organism) %>%
  summarise(Proteins = n_distinct(Entry))
print(msfraggerdia_lib_precursors)
print(msfraggerdia_lib_proteins)

lib1 <- read_tsv("msfraggerdiadda\\library.tsv", col_select = c(ProteinId, PeptideSequence, ModifiedPeptideSequence, PrecursorCharge))
lib2 <- lib1 %>%
  unique() %>%
  filter(!(PeptideSequence %in% human_ecoli_overlap_peptide$V1)) %>%
  mutate(PrecursorId = paste0(ModifiedPeptideSequence, PrecursorCharge)) %>%
  separate(col = ProteinId, into = "Entry", sep = ";", remove = TRUE, extra = "drop") %>%
  inner_join(annotation)
msfraggerdiadda_lib_precursors <- lib2 %>%
  group_by(Organism) %>%
  summarise(Precursors = n_distinct(PrecursorId))
msfraggerdiadda_lib_proteins <- lib2 %>%
  group_by(Organism) %>%
  summarise(Proteins = n_distinct(Entry))
print(msfraggerdiadda_lib_precursors)
print(msfraggerdiadda_lib_proteins)

lib1 <- read_tsv("diann\\diann-output-lib.tsv", col_select = c(ProteinGroup, PeptideSequence, PGQValue, ModifiedPeptide, PrecursorCharge))
lib2 <- lib1 %>%
  filter(PGQValue < 0.01) %>%
  select(-c(PGQValue)) %>%
  unique() %>%
  filter(!(PeptideSequence %in% human_ecoli_overlap_peptide$V1)) %>%
  mutate(PrecursorId = paste0(ModifiedPeptide, PrecursorCharge)) %>%
  separate(col = ProteinGroup, into = "Entry", sep = ";", remove = TRUE, extra = "drop") %>%
  inner_join(annotation)
diann_lib_precursors <- lib2 %>%
  group_by(Organism) %>%
  summarise(Precursors = n_distinct(PrecursorId))
diann_lib_proteins <- lib2 %>%
  group_by(Organism) %>%
  summarise(Proteins = n_distinct(Entry))
print(diann_lib_precursors)
print(diann_lib_proteins)


#### data formating etc ####

# load data
msfraggerdia <- read_tsv("msfraggerdia\\diann-output.tsv", na = c("0", ""), col_select =  c(Protein.Group, Proteotypic, Q.Value, Global.Q.Value, PG.Q.Value, Global.PG.Q.Value, Run, Precursor.Id, Precursor.Quantity, Lib.Q.Value, Lib.PG.Q.Value, Stripped.Sequence))
msfraggerdia_experimental_spectra <- read_tsv("msfraggerdia_experimental_spectra\\diann-output.tsv", na = c("0", ""), col_select =  c(Protein.Group, Proteotypic, Q.Value, Global.Q.Value, PG.Q.Value, Global.PG.Q.Value, Run, Precursor.Id, Precursor.Quantity, Lib.Q.Value, Lib.PG.Q.Value, Stripped.Sequence))
msfraggerdia_hybrid <- read_tsv("msfraggerdiadda\\diann-output.tsv", na = c("0", ""), col_select =  c(Protein.Group, Proteotypic, Q.Value, Global.Q.Value, PG.Q.Value, Global.PG.Q.Value, Run, Precursor.Id, Precursor.Quantity, Lib.Q.Value, Lib.PG.Q.Value, Stripped.Sequence))
msfraggerdia_hybrid_experimental_spectra <- read_tsv("msfraggerdiadda_experimental_spectra\\diann-output.tsv", na = c("0", ""), col_select =  c(Protein.Group, Proteotypic, Q.Value, Global.Q.Value, PG.Q.Value, Global.PG.Q.Value, Run, Precursor.Id, Precursor.Quantity, Lib.Q.Value, Lib.PG.Q.Value, Stripped.Sequence))
diann <- read_tsv("diann\\diann-output.tsv", na = c("0", ""), col_select =  c(Protein.Group, Proteotypic, Q.Value, Global.Q.Value, PG.Q.Value, Global.PG.Q.Value, Run, Precursor.Id, Precursor.Quantity, Lib.Q.Value, Lib.PG.Q.Value, Stripped.Sequence))
spectronaut_14 <- read_tsv("spectronaut\\14\\20230227_181426_directDIA_LymphEcoli_Report.tsv", na = c("NaN", "NA", ""), col_select = c(R.FileName, PG.ProteinAccessions, EG.Qvalue, EG.ModifiedSequence, FG.Charge, PG.Qvalue, "PG.QValue (Run-Wise)", "EG.TotalQuantity (Settings)", PEP.StrippedSequence))
spectronaut_17 <- read_tsv("spectronaut\\17\\20230227_180032_SN17_RealDilutionSeries_Report.tsv", na = c("NaN", "NA", ""), col_select = c(R.FileName, PG.ProteinAccessions, EG.Qvalue, EG.ModifiedSequence, FG.Charge, PG.Qvalue, "PG.QValue (Run-Wise)", "EG.TotalQuantity (Settings)", PEP.StrippedSequence))


#condition annotation
conditions_unique <- data.frame(Run = c("Lymph_001", "Lymph_002", "Lymph_004", "Lymph_0100", "Lymph_016", "Lymph_017", "Lymph_029", "Lymph_030", "Lymph_031", "Lymph_043", "Lymph_044", "Lymph_045", "Lymph_056", "Lymph_057", "Lymph_058", "Lymph_060", "Lymph_075", "Lymph_076", "Lymph_077", "Lymph_078", "Lymph_079", "Lymph_098", "Lymph_130", "Lymph_Ecoli_1-12_008", "Lymph_Ecoli_1-12_009", "Lymph_Ecoli_1-12_010", "Lymph_Ecoli_1-12_022", "Lymph_Ecoli_1-12_023", "Lymph_Ecoli_1-12_024", "Lymph_Ecoli_1-12_036", "Lymph_Ecoli_1-12_037", "Lymph_Ecoli_1-12_038", "Lymph_Ecoli_1-12_050", "Lymph_Ecoli_1-12_051", "Lymph_Ecoli_1-12_052", "Lymph_Ecoli_1-12_066", "Lymph_Ecoli_1-12_067", "Lymph_Ecoli_1-12_068", "Lymph_Ecoli_1-12_069", "Lymph_Ecoli_1-12_088", "Lymph_Ecoli_1-12_089", "Lymph_Ecoli_1-12_090", "Lymph_Ecoli_1-12_091", "Lymph_Ecoli_1-12_092", "Lymph_Ecoli_1-12_103", "Lymph_Ecoli_1-12_104", "Lymph_Ecoli_1-25_005", "Lymph_Ecoli_1-25_006", "Lymph_Ecoli_1-25_007", "Lymph_Ecoli_1-25_018", "Lymph_Ecoli_1-25_019", "Lymph_Ecoli_1-25_020", "Lymph_Ecoli_1-25_032", "Lymph_Ecoli_1-25_033", "Lymph_Ecoli_1-25_034", "Lymph_Ecoli_1-25_046", "Lymph_Ecoli_1-25_047", "Lymph_Ecoli_1-25_049", "Lymph_Ecoli_1-25_061", "Lymph_Ecoli_1-25_063", "Lymph_Ecoli_1-25_064", "Lymph_Ecoli_1-25_065", "Lymph_Ecoli_1-25_080", "Lymph_Ecoli_1-25_082", "Lymph_Ecoli_1-25_083", "Lymph_Ecoli_1-25_085", "Lymph_Ecoli_1-25_087", "Lymph_Ecoli_1-25_101", "Lymph_Ecoli_1-25_102", "Lymph_Ecoli_1-6_011", "Lymph_Ecoli_1-6_012", "Lymph_Ecoli_1-6_013", "Lymph_Ecoli_1-6_026", "Lymph_Ecoli_1-6_027", "Lymph_Ecoli_1-6_028", "Lymph_Ecoli_1-6_039", "Lymph_Ecoli_1-6_041", "Lymph_Ecoli_1-6_042", "Lymph_Ecoli_1-6_053", "Lymph_Ecoli_1-6_054", "Lymph_Ecoli_1-6_055", "Lymph_Ecoli_1-6_070", "Lymph_Ecoli_1-6_071", "Lymph_Ecoli_1-6_072", "Lymph_Ecoli_1-6_074", "Lymph_Ecoli_1-6_093", "Lymph_Ecoli_1-6_094", "Lymph_Ecoli_1-6_095", "Lymph_Ecoli_1-6_096", "Lymph_Ecoli_1-6_097", "Lymph_Ecoli_1-6_105", "Lymph_Ecoli_1-6_106"))
conditions_unique$Condition <- "Lymphnode"
conditions_unique[grepl("1-6", conditions_unique$Run), grepl("Condition", colnames(conditions_unique))] <- "1-06"
conditions_unique[grepl("1-12", conditions_unique$Run), grepl("Condition", colnames(conditions_unique))] <- "1-12"
conditions_unique[grepl("1-25", conditions_unique$Run), grepl("Condition", colnames(conditions_unique))] <- "1-25"


#### precursor level ####
# concentrate on relevant columns:
msfraggerdia <- msfraggerdia %>%
  filter(!(Stripped.Sequence %in% human_ecoli_overlap_peptide$V1)) %>%
  separate(col = Protein.Group, into = "Entry", sep = ";", remove = TRUE, extra = "drop") %>%
  inner_join(annotation) %>%
  inner_join(conditions_unique)

msfraggerdia_experimental_spectra <- msfraggerdia_experimental_spectra %>%
  filter(!(Stripped.Sequence %in% human_ecoli_overlap_peptide$V1)) %>%
  separate(col = Protein.Group, into = "Entry", sep = ";", remove = TRUE, extra = "drop") %>%
  inner_join(annotation) %>%
  inner_join(conditions_unique)

msfraggerdia_hybrid <- msfraggerdia_hybrid %>%
  filter(!(Stripped.Sequence %in% human_ecoli_overlap_peptide$V1)) %>%
  separate(col = Protein.Group, into = "Entry", sep = ";", remove = TRUE, extra = "drop") %>%
  inner_join(annotation) %>%
  inner_join(conditions_unique)

msfraggerdia_hybrid_experimental_spectra <- msfraggerdia_hybrid_experimental_spectra %>%
  filter(!(Stripped.Sequence %in% human_ecoli_overlap_peptide$V1)) %>%
  separate(col = Protein.Group, into = "Entry", sep = ";", remove = TRUE, extra = "drop") %>%
  inner_join(annotation) %>%
  inner_join(conditions_unique)

diann <- diann %>%
  filter(!(Stripped.Sequence %in% human_ecoli_overlap_peptide$V1)) %>%
  separate(col = Protein.Group, into = "Entry", sep = ";", remove = TRUE, extra = "drop") %>%
  inner_join(annotation) %>%
  inner_join(conditions_unique)

spectronaut_14 <- spectronaut_14  %>%
  filter(!(PEP.StrippedSequence %in% human_ecoli_overlap_peptide$V1)) %>%
  rename(PG.QValue.Run_Wise = "PG.QValue (Run-Wise)", Precursor.Quantity = "EG.TotalQuantity (Settings)", Run = R.FileName) %>%
  mutate(Precursor.Id = paste0(EG.ModifiedSequence, FG.Charge)) %>%
  separate(col = PG.ProteinAccessions, into = "Entry", sep = ";", remove = TRUE, extra = "drop") %>%
  inner_join(annotation) %>%
  inner_join(conditions_unique)

spectronaut_17 <- spectronaut_17  %>%
  filter(!(PEP.StrippedSequence %in% human_ecoli_overlap_peptide$V1)) %>%
  rename(PG.QValue.Run_Wise = "PG.QValue (Run-Wise)", Precursor.Quantity = "EG.TotalQuantity (Settings)", Run = R.FileName) %>%
  mutate(Precursor.Id = paste0(EG.ModifiedSequence, FG.Charge)) %>%
  separate(col = PG.ProteinAccessions, into = "Entry", sep = ";", remove = TRUE, extra = "drop") %>%
  inner_join(annotation) %>%
  inner_join(conditions_unique)


# concentrate on relevant columns:
prec_msfraggerdia <- msfraggerdia %>%
  filter(Q.Value < 0.01) %>%
  filter(Global.Q.Value < 0.01) %>%
  filter(PG.Q.Value < 0.01) %>%
  filter(Global.PG.Q.Value < 0.01) %>%
  unique()

prec_msfraggerdia_experimental_spectra <- msfraggerdia_experimental_spectra %>%
  filter(Q.Value < 0.01) %>%
  filter(Global.Q.Value < 0.01) %>%
  filter(PG.Q.Value < 0.01) %>%
  filter(Global.PG.Q.Value < 0.01) %>%
  unique()

prec_msfraggerdia_hybrid <- msfraggerdia_hybrid %>%
  filter(Q.Value < 0.01) %>%
  filter(Global.Q.Value < 0.01) %>%
  filter(PG.Q.Value < 0.01) %>%
  filter(Global.PG.Q.Value < 0.01) %>%
  unique()

prec_msfraggerdia_hybrid_experimental_spectra <- msfraggerdia_hybrid_experimental_spectra %>%
  filter(Q.Value < 0.01) %>%
  filter(Global.Q.Value < 0.01) %>%
  filter(PG.Q.Value < 0.01) %>%
  filter(Global.PG.Q.Value < 0.01) %>%
  unique()

prec_diann <- diann %>%
  filter(Q.Value < 0.01) %>%
  filter(Global.Q.Value < 0.01) %>%
  filter(PG.Q.Value < 0.01) %>%
  filter(Global.PG.Q.Value < 0.01) %>%
  unique()

prec_spectronaut_14 <- spectronaut_14 %>%
  filter(PG.Qvalue < 0.01) %>%
  # filter(PG.QValue.Run_Wise < 0.01) %>%  # Spectronaut 14 does not have run-wise q-value
  filter(EG.Qvalue < 0.01) %>%
  unique()

prec_spectronaut_17 <- spectronaut_17 %>%
  filter(PG.Qvalue < 0.01) %>%
  filter(PG.QValue.Run_Wise < 0.01) %>%
  filter(EG.Qvalue < 0.01) %>%
  unique()


# one big list
prec_list <- list()
prec_list$`Spectronaut 14` <- prec_spectronaut_14
prec_list$`Spectronaut 17` <- prec_spectronaut_17
prec_list$`DIA-NN lib-free` <- prec_diann
prec_list$`FP-MSF` <- prec_msfraggerdia
prec_list$`FP-MSF hybrid` <- prec_msfraggerdia_hybrid


#combine to one dataframe
prec_df <- bind_rows(prec_list, .id = "Origin") %>%
  drop_na(Organism) %>%
  drop_na(Condition)

prec_df$Condition <- factor(prec_df$Condition)


# peptide IDs per software solution separated by conditions
prec_IDs <- prec_df %>%
  group_by(Origin, Condition, Run, Organism) %>%
  summarise(Precursors = n_distinct(Precursor.Id))

prec_IDs$Origin <- factor(prec_IDs$Origin)

ggplot(prec_IDs, aes(x = Precursors, y = Origin, color = Condition)) +
  geom_boxplot(size = 0.5) +
  facet_grid(~ Organism, scales = "free") +
  theme(text = element_text(size = 20), axis.title.y = element_blank()) +
  scale_y_discrete(limits = c("FP-MSF hybrid", "FP-MSF", "DIA-NN lib-free", "Spectronaut 17", "Spectronaut 14")) +
  geom_vline(xintercept = 0, size=0.2, colour="red")

ggsave("benchmark_boxplot_precursors.pdf", width = 15, height = 6, units = "in")


## print numbers
tool <- "Spectronaut 14"
yy <- prec_spectronaut_14[prec_spectronaut_14$Condition == "Lymphnode", ] %>%
  group_by(Organism) %>%
  summarise(Precursors = n_distinct(Precursor.Id))
a <- yy[yy$Organism == "Escherichia coli (strain K12)", ]$Precursors
b <- yy[yy$Organism == "Homo sapiens", ]$Precursors
sprintf("%s, %.1f, %.1f, %.3f, NA", tool, a, b, a * 100 / b)

tool <- "Spectronaut 17"
yy <- prec_spectronaut_17[prec_spectronaut_17$Condition == "Lymphnode", ] %>%
  group_by(Organism) %>%
  summarise(Precursors = n_distinct(Precursor.Id))
a <- yy[yy$Organism == "Escherichia coli (strain K12)", ]$Precursors
b <- yy[yy$Organism == "Homo sapiens", ]$Precursors
sprintf("%s, %.1f, %.1f, %.3f, NA", tool, a, b, a * 100 / b)

tool <- "DIA-NN lib-free"
yy <- prec_diann[prec_diann$Condition == "Lymphnode", ] %>%
  group_by(Organism) %>%
  summarise(Precursors = n_distinct(Precursor.Id))
a <- yy[yy$Organism == "Escherichia coli (strain K12)", ]$Precursors
b <- yy[yy$Organism == "Homo sapiens", ]$Precursors
c <- diann_lib_precursors[diann_lib_precursors$Organism == "Escherichia coli (strain K12)", ]$Precursors
d <- diann_lib_precursors[diann_lib_precursors$Organism == "Homo sapiens", ]$Precursors
sprintf("%s, %.1f, %.1f, %.3f, %.3f", tool, a, b, a * 100 / b, a * (d - b) * 100 / (b * c))

tool <- "FP-MSF"
yy <- prec_msfraggerdia[prec_msfraggerdia$Condition == "Lymphnode", ] %>%
  group_by(Organism) %>%
  summarise(Precursors = n_distinct(Precursor.Id))
a <- yy[yy$Organism == "Escherichia coli (strain K12)", ]$Precursors
b <- yy[yy$Organism == "Homo sapiens", ]$Precursors
c <- msfraggerdia_lib_precursors[msfraggerdia_lib_precursors$Organism == "Escherichia coli (strain K12)", ]$Precursors
d <- msfraggerdia_lib_precursors[msfraggerdia_lib_precursors$Organism == "Homo sapiens", ]$Precursors
sprintf("%s, %.1f, %.1f, %.3f, %.3f", tool, a, b, a * 100 / b, a * (d - b) * 100 / (b * c))

tool <- "FP-MSF experimental spectra"
yy <- prec_msfraggerdia_experimental_spectra[prec_msfraggerdia_experimental_spectra$Condition == "Lymphnode", ] %>%
  group_by(Organism) %>%
  summarise(Precursors = n_distinct(Precursor.Id))
a <- yy[yy$Organism == "Escherichia coli (strain K12)", ]$Precursors
b <- yy[yy$Organism == "Homo sapiens", ]$Precursors
c <- msfraggerdia_lib_precursors[msfraggerdia_lib_precursors$Organism == "Escherichia coli (strain K12)", ]$Precursors
d <- msfraggerdia_lib_precursors[msfraggerdia_lib_precursors$Organism == "Homo sapiens", ]$Precursors
sprintf("%s, %.1f, %.1f, %.3f, %.3f", tool, a, b, a * 100 / b, a * (d - b) * 100 / (b * c))

tool <- "FP-MSF hybrid"
yy <- prec_msfraggerdia_hybrid[prec_msfraggerdia_hybrid$Condition == "Lymphnode", ] %>%
  group_by(Organism) %>%
  summarise(Precursors = n_distinct(Precursor.Id))
a <- yy[yy$Organism == "Escherichia coli (strain K12)", ]$Precursors
b <- yy[yy$Organism == "Homo sapiens", ]$Precursors
c <- msfraggerdiadda_lib_precursors[msfraggerdiadda_lib_precursors$Organism == "Escherichia coli (strain K12)", ]$Precursors
d <- msfraggerdiadda_lib_precursors[msfraggerdiadda_lib_precursors$Organism == "Homo sapiens", ]$Precursors
sprintf("%s, %.1f, %.1f, %.3f, %.3f", tool, a, b, a * 100 / b, a * (d - b) * 100 / (b * c))

tool <- "FP-MSF hybrid experimental spectra"
yy <- prec_msfraggerdia_hybrid_experimental_spectra[prec_msfraggerdia_hybrid_experimental_spectra$Condition == "Lymphnode", ] %>%
  group_by(Organism) %>%
  summarise(Precursors = n_distinct(Precursor.Id))
a <- yy[yy$Organism == "Escherichia coli (strain K12)", ]$Precursors
b <- yy[yy$Organism == "Homo sapiens", ]$Precursors
c <- msfraggerdiadda_lib_precursors[msfraggerdiadda_lib_precursors$Organism == "Escherichia coli (strain K12)", ]$Precursors
d <- msfraggerdiadda_lib_precursors[msfraggerdiadda_lib_precursors$Organism == "Homo sapiens", ]$Precursors
sprintf("%s, %.1f, %.1f, %.3f, %.3f", tool, a, b, a * 100 / b, a * (d - b) * 100 / (b * c))
