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

rm(list = ls())
library(LFQbench)
library(readr)
library(dplyr)
library(stringr)
library(tidyr)


my_process_1 <- function(input_path, tool, annotation2_path, is_spectronaut) {
  if (is_spectronaut) {
    annotation2 <- read_tsv(annotation2_path, col_select = c("PG.ProteinGroups", "PEP.StrippedSequence", "EG.ModifiedSequence", "FG.Charge")) %>%
      rename(Protein.Group = PG.ProteinGroups, Stripped.Sequence = PEP.StrippedSequence) %>%
      mutate(Precursor.Id = paste0(EG.ModifiedSequence, FG.Charge))
  } else {
    annotation2 <- read_tsv(annotation2_path, col_select = c("Protein.Group", "Stripped.Sequence", "Precursor.Id"))
  }

  annotation2 <- annotation2 %>%
    distinct(Precursor.Id, .keep_all = TRUE) %>%
    separate(col = Protein.Group, into = "Entry", sep = ";", remove = TRUE, extra = "drop") %>%
    inner_join(annotation)

  input_table <- read_tsv(file = input_path)

  colnames(input_table)[2:length(colnames(input_table))] <- gsub("$", " intensity", colnames(input_table)[2:length(colnames(input_table))])

  input_table <- input_table %>%
    rename(Precursor.Id = `...1`) %>%
    select(c("Precursor.Id", contains("1-6"), contains("1-12"))) %>%
    inner_join(annotation2) %>%
    filter(!(Stripped.Sequence %in% human_ecoli_overlap_peptide$V1)) %>%
    mutate(Protein = paste0("sp|", Precursor.Id, "|", `Entry name`))

  colnames(input_table) <- sub("-", "_", colnames(input_table))
  write_tsv(input_table, str_c(srcDir, str_c(tool, ".tsv")))
}


setwd("G:\\Dropbox\\papers_Fengchao\\msfragger_dia\\script\\results\\benchmark\\")

annotation <- read_tsv("G:\\Dropbox\\papers_Fengchao\\msfragger_dia\\script\\code\\MSFragger-DIA-manuscript\\uniprot.tab", col_select = c("Entry", "Entry name"))

human_ecoli_overlap_peptide <- read.delim("human_ecoli_overlap_peptides.txt", header = FALSE, sep = "\t")

srcDir <- "./lfqbench_benchmark_2to1_precursor/"
unlink(srcDir, recursive = TRUE)
mkdir(srcDir)

my_process_1("spectronaut\\14\\precursor_maxlfq.tsv", "my_spectronaut_14", "spectronaut\\14\\20230227_181426_directDIA_LymphEcoli_Report.tsv", TRUE)
my_process_1("spectronaut\\17\\precursor_maxlfq.tsv", "my_spectronaut_17", "spectronaut\\17\\20230227_180032_SN17_RealDilutionSeries_Report.tsv", TRUE)
my_process_1("diann\\precursor_maxlfq.tsv", "diann", "diann\\diann-output.tsv", FALSE)
my_process_1("msfraggerdia\\precursor_maxlfq.tsv", "msfraggerdia", "msfraggerdia\\diann-output.tsv", FALSE)
my_process_1("msfraggerdiadda\\precursor_maxlfq.tsv", "msfragger_hybrid", "msfraggerdiadda\\diann-output.tsv", FALSE)


sampleComposition <- data.frame(
  species = c("HUMAN", "ECOLI"),
  A = c(1, 2),
  B = c(1, 1)
)

dataSets <- data.frame(
  "HYE124" = c("Lymph_Ecoli_1_6_011",
               "Lymph_Ecoli_1_6_012",
               "Lymph_Ecoli_1_6_013",
               "Lymph_Ecoli_1_6_026",
               "Lymph_Ecoli_1_6_027",
               "Lymph_Ecoli_1_6_028",
               "Lymph_Ecoli_1_6_039",
               "Lymph_Ecoli_1_6_041",
               "Lymph_Ecoli_1_6_042",
               "Lymph_Ecoli_1_6_053",
               "Lymph_Ecoli_1_6_054",
               "Lymph_Ecoli_1_6_055",
               "Lymph_Ecoli_1_6_070",
               "Lymph_Ecoli_1_6_071",
               "Lymph_Ecoli_1_6_072",
               "Lymph_Ecoli_1_6_074",
               "Lymph_Ecoli_1_6_093",
               "Lymph_Ecoli_1_6_094",
               "Lymph_Ecoli_1_6_095",
               "Lymph_Ecoli_1_6_096",
               "Lymph_Ecoli_1_6_097",
               "Lymph_Ecoli_1_6_105",
               "Lymph_Ecoli_1_6_106",
               "Lymph_Ecoli_1_12_008",
               "Lymph_Ecoli_1_12_009",
               "Lymph_Ecoli_1_12_010",
               "Lymph_Ecoli_1_12_022",
               "Lymph_Ecoli_1_12_023",
               "Lymph_Ecoli_1_12_024",
               "Lymph_Ecoli_1_12_036",
               "Lymph_Ecoli_1_12_037",
               "Lymph_Ecoli_1_12_038",
               "Lymph_Ecoli_1_12_050",
               "Lymph_Ecoli_1_12_051",
               "Lymph_Ecoli_1_12_052",
               "Lymph_Ecoli_1_12_066",
               "Lymph_Ecoli_1_12_067",
               "Lymph_Ecoli_1_12_068",
               "Lymph_Ecoli_1_12_069",
               "Lymph_Ecoli_1_12_088",
               "Lymph_Ecoli_1_12_089",
               "Lymph_Ecoli_1_12_090",
               "Lymph_Ecoli_1_12_091",
               "Lymph_Ecoli_1_12_092",
               "Lymph_Ecoli_1_12_103",
               "Lymph_Ecoli_1_12_104"
  ),
  row.names = c("A1", "A2", "A3", "A4", "A5", "A6", "A7", "A8", "A9", "A10", "A11", "A12", "A13", "A14", "A15", "A16", "A17", "A18", "A19", "A20", "A21", "A22", "A23", "B1", "B2", "B3", "B4", "B5", "B6", "B7", "B8", "B9", "B10", "B11", "B12", "B13", "B14", "B15", "B16", "B17", "B18", "B19", "B20", "B21", "B22", "B23")
)

speciesTags <- list(
  HUMAN = "_HUMAN",
  ECOLI = "_ECOLI"
)

LFQbench.initConfiguration(
  SampleComposition = sampleComposition
)
FSWE.initConfiguration(
  injectionNames = dataSets,
  speciesTags = speciesTags
)

LFQbench.setDataRootFolder(
  rootFolder = srcDir,
  createSubfolders = FALSE
)

FSWE.addSoftwareConfiguration(
  softwareName = "my_spectronaut",
  input_format = "wide",
  input.extension = "*.tsv$",
  nastrings = "NA",
  protein_input = TRUE,
  quantitative.var = "intensity",
  quantitative.var.tag = "intensity",
  protein.var = "Protein",
)

FSWE.addSoftwareConfiguration(
  softwareName = "diann",
  input_format = "wide",
  input.extension = "*.tsv$",
  nastrings = "NA",
  protein_input = TRUE,
  quantitative.var = "intensity",
  quantitative.var.tag = "intensity",
  protein.var = "Protein",
)

FSWE.addSoftwareConfiguration(
  softwareName = "msfraggerdia",
  input_format = "wide",
  input.extension = "*.tsv$",
  nastrings = "NA",
  protein_input = TRUE,
  quantitative.var = "intensity",
  quantitative.var.tag = "intensity",
  protein.var = "Protein",
)

FSWE.addSoftwareConfiguration(
  softwareName = "msfragger_hybrid",
  input_format = "wide",
  input.extension = "*.tsv$",
  nastrings = "NA",
  protein_input = TRUE,
  quantitative.var = "intensity",
  quantitative.var.tag = "intensity",
  protein.var = "Protein",
)

inputFiles <- list.files(
  path = LFQbench.Config$DataRootFolder,
  pattern = ".+\\.tsv"
)

nix <- sapply(
  inputFiles,
  FSWE.generateReports,
  softwareSource = "guess",
  keep_original_names = T,
  singleHits = F,
  plotHistogram = T,
  plotHistNAs = T,
  reportSequences = F
)

# some configuration changes for beautifying plots
LFQbench.changeConfiguration(
  PlotPointSize = 0.5,
  LogRatioPlotRange = c(-5, 5),
  LogRatioValidityRangeSDFactor = 50
)
# run batch analysis and keep result set
res <- LFQbench.batchProcessRootFolder()
