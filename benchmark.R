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

setwd("G:\\Dropbox\\papers_Fengchao\\msfragger_dia\\script\\code\\MSFragger-DIA-manuscript\\")


# load annotation data
annotation <- read_tsv("G:\\Dropbox\\papers_Fengchao\\msfragger_dia\\script\\code\\MSFragger-DIA-manuscript\\uniprot.tab", col_select = c("Entry", "Organism"))


#### list library stats ####
lib1 <- read_tsv("G:\\Dropbox\\papers_Fengchao\\msfragger_dia\\script\\results\\benchmark\\msfraggerdia\\library.tsv", col_select = c(ProteinId, ModifiedPeptideSequence, PrecursorCharge))
lib2 <- lib1 %>%
  unique() %>%
  mutate(PrecursorId = paste0(ModifiedPeptideSequence, PrecursorCharge)) %>%
  separate(col = ProteinId, into = "Entry", sep = ";") %>%
  inner_join(annotation)
msfraggerdia_lib_precursors <- lib2 %>%
  group_by(Organism) %>%
  summarise(Precursors = n_distinct(PrecursorId))
msfraggerdia_lib_proteins <- lib2 %>%
  group_by(Organism) %>%
  summarise(Proteins = n_distinct(Entry))

lib1 <- read_tsv("G:\\Dropbox\\papers_Fengchao\\msfragger_dia\\script\\results\\benchmark\\msfraggerdiadda\\library.tsv", col_select = c(ProteinId, ModifiedPeptideSequence, PrecursorCharge))
lib2 <- lib1 %>%
  unique() %>%
  mutate(PrecursorId = paste0(ModifiedPeptideSequence, PrecursorCharge)) %>%
  separate(col = ProteinId, into = "Entry", sep = ";") %>%
  inner_join(annotation)
msfraggerdiadda_lib_precursors <- lib2 %>%
  group_by(Organism) %>%
  summarise(Precursors = n_distinct(PrecursorId))
msfraggerdiadda_lib_proteins <- lib2 %>%
  group_by(Organism) %>%
  summarise(Proteins = n_distinct(Entry))

lib1 <- read_tsv("G:\\Dropbox\\papers_Fengchao\\msfragger_dia\\script\\results\\benchmark\\diann\\diann-output-lib.tsv", col_select = c(ProteinGroup, PGQValue, ModifiedPeptide, PrecursorCharge))
lib2 <- lib1 %>%
  filter(PGQValue < 0.01) %>%
  unique() %>%
  mutate(PrecursorId = paste0(ModifiedPeptide, PrecursorCharge)) %>%
  separate(col = ProteinGroup, into = "Entry", sep = ";") %>%
  inner_join(annotation)
diann_lib_precursors <- lib2 %>%
  group_by(Organism) %>%
  summarise(Precursors = n_distinct(PrecursorId))
diann_lib_proteins <- lib2 %>%
  group_by(Organism) %>%
  summarise(Proteins = n_distinct(Entry))


#### data formating etc ####

# load data
msfraggerdia <- read_tsv("G:\\Dropbox\\papers_Fengchao\\msfragger_dia\\script\\results\\benchmark\\msfraggerdia\\diann-output.tsv", na = c("NA", "0", ""), col_select =  c(Protein.Group, Proteotypic, Q.Value, Global.Q.Value, PG.Q.Value, Global.PG.Q.Value, Run, Precursor.Id, Precursor.Quantity, Lib.Q.Value, Lib.PG.Q.Value))
msfraggerdia_experimental_spectra <- read_tsv("G:\\Dropbox\\papers_Fengchao\\msfragger_dia\\script\\results\\benchmark\\msfraggerdia_experimental_spectra\\diann-output.tsv", na = c("NA", "0", ""), col_select =  c(Protein.Group, Proteotypic, Q.Value, Global.Q.Value, PG.Q.Value, Global.PG.Q.Value, Run, Precursor.Id, Precursor.Quantity, Lib.Q.Value, Lib.PG.Q.Value))
msfraggerdia_hybrid <- read_tsv("G:\\Dropbox\\papers_Fengchao\\msfragger_dia\\script\\results\\benchmark\\msfraggerdiadda\\diann-output.tsv", na = c("NA", "0", ""), col_select =  c(Protein.Group, Proteotypic, Q.Value, Global.Q.Value, PG.Q.Value, Global.PG.Q.Value, Run, Precursor.Id, Precursor.Quantity, Lib.Q.Value, Lib.PG.Q.Value))
msfraggerdia_hybrid_experimental_spectra <- read_tsv("G:\\Dropbox\\papers_Fengchao\\msfragger_dia\\script\\results\\benchmark\\msfraggerdiadda_experimental_spectra\\diann-output.tsv", na = c("NA", "0", ""), col_select =  c(Protein.Group, Proteotypic, Q.Value, Global.Q.Value, PG.Q.Value, Global.PG.Q.Value, Run, Precursor.Id, Precursor.Quantity, Lib.Q.Value, Lib.PG.Q.Value))
diann <- read_tsv("G:\\Dropbox\\papers_Fengchao\\msfragger_dia\\script\\results\\benchmark\\diann\\diann-output.tsv", na = c("NA", "0", ""), col_select =  c(Protein.Group, Proteotypic, Q.Value, Global.Q.Value, PG.Q.Value, Global.PG.Q.Value, Run, Precursor.Id, Precursor.Quantity, Lib.Q.Value, Lib.PG.Q.Value))
Spectronaut_qVal <- read_tsv("G:\\Dropbox\\papers_Fengchao\\msfragger_dia\\script\\results\\benchmark\\spectronaut\\Spectronaut_DirectDIAV15_qVal.xls", na = c("NA", "0", ""), col_select = c(R.FileName, PEP.IsProteotypic, PG.ProteinAccessions, EG.Qvalue, EG.ModifiedSequence, FG.Charge, PG.Qvalue, "PG.QValue (Run-Wise)", "EG.TotalQuantity (Settings)"))


#condition annotation
conditions_unique <- data.frame(Run = c("Lymph_001", "Lymph_002", "Lymph_004", "Lymph_0100", "Lymph_016", "Lymph_017", "Lymph_029", "Lymph_030", "Lymph_031", "Lymph_043", "Lymph_044", "Lymph_045", "Lymph_056", "Lymph_057", "Lymph_058", "Lymph_060", "Lymph_075", "Lymph_076", "Lymph_077", "Lymph_078", "Lymph_079", "Lymph_098", "Lymph_130", "Lymph_Ecoli_1-12_008", "Lymph_Ecoli_1-12_009", "Lymph_Ecoli_1-12_010", "Lymph_Ecoli_1-12_022", "Lymph_Ecoli_1-12_023", "Lymph_Ecoli_1-12_024", "Lymph_Ecoli_1-12_036", "Lymph_Ecoli_1-12_037", "Lymph_Ecoli_1-12_038", "Lymph_Ecoli_1-12_050", "Lymph_Ecoli_1-12_051", "Lymph_Ecoli_1-12_052", "Lymph_Ecoli_1-12_066", "Lymph_Ecoli_1-12_067", "Lymph_Ecoli_1-12_068", "Lymph_Ecoli_1-12_069", "Lymph_Ecoli_1-12_088", "Lymph_Ecoli_1-12_089", "Lymph_Ecoli_1-12_090", "Lymph_Ecoli_1-12_091", "Lymph_Ecoli_1-12_092", "Lymph_Ecoli_1-12_103", "Lymph_Ecoli_1-12_104", "Lymph_Ecoli_1-25_005", "Lymph_Ecoli_1-25_006", "Lymph_Ecoli_1-25_007", "Lymph_Ecoli_1-25_018", "Lymph_Ecoli_1-25_019", "Lymph_Ecoli_1-25_020", "Lymph_Ecoli_1-25_032", "Lymph_Ecoli_1-25_033", "Lymph_Ecoli_1-25_034", "Lymph_Ecoli_1-25_046", "Lymph_Ecoli_1-25_047", "Lymph_Ecoli_1-25_049", "Lymph_Ecoli_1-25_061", "Lymph_Ecoli_1-25_063", "Lymph_Ecoli_1-25_064", "Lymph_Ecoli_1-25_065", "Lymph_Ecoli_1-25_080", "Lymph_Ecoli_1-25_082", "Lymph_Ecoli_1-25_083", "Lymph_Ecoli_1-25_085", "Lymph_Ecoli_1-25_087", "Lymph_Ecoli_1-25_101", "Lymph_Ecoli_1-25_102", "Lymph_Ecoli_1-6_011", "Lymph_Ecoli_1-6_012", "Lymph_Ecoli_1-6_013", "Lymph_Ecoli_1-6_026", "Lymph_Ecoli_1-6_027", "Lymph_Ecoli_1-6_028", "Lymph_Ecoli_1-6_039", "Lymph_Ecoli_1-6_041", "Lymph_Ecoli_1-6_042", "Lymph_Ecoli_1-6_053", "Lymph_Ecoli_1-6_054", "Lymph_Ecoli_1-6_055", "Lymph_Ecoli_1-6_070", "Lymph_Ecoli_1-6_071", "Lymph_Ecoli_1-6_072", "Lymph_Ecoli_1-6_074", "Lymph_Ecoli_1-6_093", "Lymph_Ecoli_1-6_094", "Lymph_Ecoli_1-6_095", "Lymph_Ecoli_1-6_096", "Lymph_Ecoli_1-6_097", "Lymph_Ecoli_1-6_105", "Lymph_Ecoli_1-6_106"))
conditions_unique$Condition <- "Lymphnode"
conditions_unique[grepl("1-6", conditions_unique$Run), grepl("Condition", colnames(conditions_unique))] <- "1-06"
conditions_unique[grepl("1-12", conditions_unique$Run), grepl("Condition", colnames(conditions_unique))] <- "1-12"
conditions_unique[grepl("1-25", conditions_unique$Run), grepl("Condition", colnames(conditions_unique))] <- "1-25"


#### precursor level ####
# concentrate on relevant columns:
msfraggerdia <- msfraggerdia %>%
  separate(col = Protein.Group, into = "Entry", sep = ";") %>%
  inner_join(annotation) %>%
  inner_join(conditions_unique)

msfraggerdia_experimental_spectra <- msfraggerdia_experimental_spectra %>%
  separate(col = Protein.Group, into = "Entry", sep = ";") %>%
  inner_join(annotation) %>%
  inner_join(conditions_unique)

msfraggerdia_hybrid <- msfraggerdia_hybrid %>%
  separate(col = Protein.Group, into = "Entry", sep = ";") %>%
  inner_join(annotation) %>%
  inner_join(conditions_unique)

msfraggerdia_hybrid_experimental_spectra <- msfraggerdia_hybrid_experimental_spectra %>%
  separate(col = Protein.Group, into = "Entry", sep = ";") %>%
  inner_join(annotation) %>%
  inner_join(conditions_unique)

diann <- diann %>%
  separate(col = Protein.Group, into = "Entry", sep = ";") %>%
  inner_join(annotation) %>%
  inner_join(conditions_unique)

Spectronaut_qVal <- Spectronaut_qVal  %>%
  rename(PG.QValue.Run_Wise = "PG.QValue (Run-Wise)", Precursor.Quantity = "EG.TotalQuantity (Settings)", Run = R.FileName) %>%
  mutate(Precursor.Id = paste0(EG.ModifiedSequence, FG.Charge)) %>%
  separate(col = PG.ProteinAccessions, into = "Entry", sep = ";") %>%
  inner_join(annotation) %>%
  inner_join(conditions_unique)


# concentrate on relevant columns:
prec_msfraggerdia <- msfraggerdia %>%
  # filter(is.finite(Proteotypic)) %>%
  filter(Q.Value < 0.01) %>%
  filter(Global.Q.Value < 0.01) %>%
  filter(PG.Q.Value < 0.01) %>%
  filter(Global.PG.Q.Value < 0.01) %>%
  unique()

prec_msfraggerdia_experimental_spectra <- msfraggerdia_experimental_spectra %>%
  # filter(is.finite(Proteotypic)) %>%
  filter(Q.Value < 0.01) %>%
  filter(Global.Q.Value < 0.01) %>%
  filter(PG.Q.Value < 0.01) %>%
  filter(Global.PG.Q.Value < 0.01) %>%
  unique()

prec_msfraggerdia_hybrid <- msfraggerdia_hybrid %>%
  # filter(is.finite(Proteotypic)) %>%
  filter(Q.Value < 0.01) %>%
  filter(Global.Q.Value < 0.01) %>%
  filter(PG.Q.Value < 0.01) %>%
  filter(Global.PG.Q.Value < 0.01) %>%
  unique()

prec_msfraggerdia_hybrid_experimental_spectra <- msfraggerdia_hybrid_experimental_spectra %>%
  # filter(is.finite(Proteotypic)) %>%
  filter(Q.Value < 0.01) %>%
  filter(Global.Q.Value < 0.01) %>%
  filter(PG.Q.Value < 0.01) %>%
  filter(Global.PG.Q.Value < 0.01) %>%
  unique()

prec_diann <- diann %>%
  # filter(is.finite(Proteotypic)) %>%
  filter(Q.Value < 0.01) %>%
  filter(Global.Q.Value < 0.01) %>%
  filter(PG.Q.Value < 0.01) %>%
  filter(Global.PG.Q.Value < 0.01) %>%
  unique()

prec_Spectronaut_DirectDIA_qVal <- Spectronaut_qVal  %>%
  # filter(PEP.IsProteotypic == TRUE) %>%
  filter(PG.Qvalue < 0.01) %>%
  filter(PG.QValue.Run_Wise < 0.01) %>%
  filter(EG.Qvalue < 0.01) %>%
  unique()


# one big list
prec_list <- list()
prec_list$Spectronaut <- prec_Spectronaut_DirectDIA_qVal
prec_list$`DIA-NN` <- prec_diann
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

ggplot(prec_IDs, aes(x = Precursors , y = Origin, color = Condition))+
  geom_boxplot(size = 0.5) +
  facet_grid(~ Organism, scales = "free") +
  theme(text = element_text(size = 20), axis.title.y = element_blank()) +
  scale_y_discrete(limits = c("FP-MSF hybrid", "FP-MSF", "DIA-NN", "Spectronaut"))

ggsave("benchmark_precursors.pdf", width = 15, height = 6, units = "in")


# # precursor variance per software solution separated by condition
# prec_var <- prec_df %>%
#   group_by(Origin, Condition, Entry, Organism) %>%
#   summarise(Variance = sd(log2(Precursor.Quantity), na.rm = TRUE))
#
# ggplot(prec_var, aes(x = Variance , y = Origin, color = Condition))+
#   geom_boxplot(size = 0.5) +
#   facet_grid(~ Organism, scales = "free") +
#   theme(text = element_text(size = 20), axis.title.y = element_blank())
#
# ggsave("lumph_variance.pdf", width = 15, height = 6, units = "in")


## print numbers

# one big list
prec_list <- list()
prec_list$Spectronaut <- prec_Spectronaut_DirectDIA_qVal
prec_list$`DIA-NN` <- prec_diann
prec_list$`MSFragger` <- prec_msfraggerdia
prec_list$`MSFragger experimental spectra` <- prec_msfraggerdia_experimental_spectra
prec_list$`MSFragger hybrid` <- prec_msfraggerdia_hybrid
prec_list$`MSFragger hybrid experimental spectra` <- prec_msfraggerdia_hybrid_experimental_spectra


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


tool <- "Spectronaut"
yy <- prec_IDs[prec_IDs$Condition == "Lymphnode" & prec_IDs$Origin == tool, ]
a <- median(yy[yy$Organism == "Escherichia coli (strain K12)", ]$Precursors)
b <- median(yy[yy$Organism == "Homo sapiens", ]$Precursors)
sprintf("%s, %d, %d, %.3f, NA", tool, a, b, a * 100 / b)

tool <- "DIA-NN"
yy <- prec_IDs[prec_IDs$Condition == "Lymphnode" & prec_IDs$Origin == tool, ]
a <- median(yy[yy$Organism == "Escherichia coli (strain K12)", ]$Precursors)
b <- median(yy[yy$Organism == "Homo sapiens", ]$Precursors)
c <- diann_lib_precursors[diann_lib_precursors$Organism == "Escherichia coli (strain K12)", ]$Precursors
d <- diann_lib_precursors[diann_lib_precursors$Organism == "Homo sapiens", ]$Precursors
sprintf("%s, %d, %d, %.3f, %.3f", tool, a, b, a * 100 / b, a * (d - b) * 100 / (b * c))

tool <- "MSFragger"
yy <- prec_IDs[prec_IDs$Condition == "Lymphnode" & prec_IDs$Origin == tool, ]
a <- median(yy[yy$Organism == "Escherichia coli (strain K12)", ]$Precursors)
b <- median(yy[yy$Organism == "Homo sapiens", ]$Precursors)
c <- msfraggerdia_lib_precursors[msfraggerdia_lib_precursors$Organism == "Escherichia coli (strain K12)", ]$Precursors
d <- msfraggerdia_lib_precursors[msfraggerdia_lib_precursors$Organism == "Homo sapiens", ]$Precursors
sprintf("%s, %d, %d, %.3f, %.3f", tool, a, b, a * 100 / b, a * (d - b) * 100 / (b * c))

tool <- "MSFragger experimental spectra"
yy <- prec_IDs[prec_IDs$Condition == "Lymphnode" & prec_IDs$Origin == tool, ]
a <- median(yy[yy$Organism == "Escherichia coli (strain K12)", ]$Precursors)
b <- median(yy[yy$Organism == "Homo sapiens", ]$Precursors)
c <- msfraggerdia_lib_precursors[msfraggerdia_lib_precursors$Organism == "Escherichia coli (strain K12)", ]$Precursors
d <- msfraggerdia_lib_precursors[msfraggerdia_lib_precursors$Organism == "Homo sapiens", ]$Precursors
sprintf("%s, %d, %d, %.3f, %.3f", tool, a, b, a * 100 / b, a * (d - b) * 100 / (b * c))

tool <- "MSFragger hybrid"
yy <- prec_IDs[prec_IDs$Condition == "Lymphnode" & prec_IDs$Origin == tool, ]
a <- median(yy[yy$Organism == "Escherichia coli (strain K12)", ]$Precursors)
b <- median(yy[yy$Organism == "Homo sapiens", ]$Precursors)
c <- msfraggerdiadda_lib_precursors[msfraggerdiadda_lib_precursors$Organism == "Escherichia coli (strain K12)", ]$Precursors
d <- msfraggerdiadda_lib_precursors[msfraggerdiadda_lib_precursors$Organism == "Homo sapiens", ]$Precursors
sprintf("%s, %d, %d, %.3f, %.3f", tool, a, b, a * 100 / b, a * (d - b) * 100 / (b * c))

tool <- "MSFragger hybrid experimental spectra"
yy <- prec_IDs[prec_IDs$Condition == "Lymphnode" & prec_IDs$Origin == tool, ]
a <- median(yy[yy$Organism == "Escherichia coli (strain K12)", ]$Precursors)
b <- median(yy[yy$Organism == "Homo sapiens", ]$Precursors)
c <- msfraggerdiadda_lib_precursors[msfraggerdiadda_lib_precursors$Organism == "Escherichia coli (strain K12)", ]$Precursors
d <- msfraggerdiadda_lib_precursors[msfraggerdiadda_lib_precursors$Organism == "Homo sapiens", ]$Precursors
sprintf("%s, %d, %d, %.3f, %.3f", tool, a, b, a * 100 / b, a * (d - b) * 100 / (b * c))
