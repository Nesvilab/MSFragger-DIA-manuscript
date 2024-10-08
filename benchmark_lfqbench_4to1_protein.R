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

my_process_1 <- function(input_path, tool) {
  input_table <- read_tsv(file = input_path)

  colnames(input_table)[2:length(colnames(input_table))] <- gsub("$", " intensity", colnames(input_table)[2:length(colnames(input_table))])

  input_table <- input_table %>%
    select(c("...1", contains("1-6"), contains("1-25"))) %>%
    separate(col = ...1, into = "Entry", sep = ";", remove = TRUE, extra = "drop") %>%
    inner_join(annotation) %>%
    mutate(Protein = paste0("sp|", Entry, "|", `Entry name`))

  colnames(input_table) <- sub("-", "_", colnames(input_table))
  write_tsv(input_table, str_c(srcDir, str_c(tool, ".tsv")))
}


setwd("G:\\Dropbox\\papers_Fengchao\\msfragger_dia\\script\\results\\benchmark\\")

annotation <- read_tsv("G:\\Dropbox\\papers_Fengchao\\msfragger_dia\\script\\code\\MSFragger-DIA-manuscript\\uniprot.tab", col_select = c("Entry", "Entry name"))

srcDir <- "./lfqbench_benchmark_4to1_protein/"
unlink(srcDir, recursive = TRUE)
mkdir(srcDir)

my_process_1("spectronaut\\14\\protein_maxlfq.tsv", "my_spectronaut_14")
my_process_1("spectronaut\\17\\protein_maxlfq.tsv", "my_spectronaut_17")
my_process_1("diann\\protein_maxlfq.tsv", "diann")
my_process_1("msfraggerdia\\protein_maxlfq.tsv", "msfraggerdia")
my_process_1("msfraggerdiadda\\protein_maxlfq.tsv", "msfragger_hybrid")


sampleComposition <- data.frame(
  species = c("HUMAN", "ECOLI"),
  A = c(1, 25),
  B = c(1, 6)
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
               "Lymph_Ecoli_1_25_005",
               "Lymph_Ecoli_1_25_006",
               "Lymph_Ecoli_1_25_007",
               "Lymph_Ecoli_1_25_018",
               "Lymph_Ecoli_1_25_019",
               "Lymph_Ecoli_1_25_020",
               "Lymph_Ecoli_1_25_032",
               "Lymph_Ecoli_1_25_033",
               "Lymph_Ecoli_1_25_034",
               "Lymph_Ecoli_1_25_046",
               "Lymph_Ecoli_1_25_047",
               "Lymph_Ecoli_1_25_049",
               "Lymph_Ecoli_1_25_061",
               "Lymph_Ecoli_1_25_063",
               "Lymph_Ecoli_1_25_064",
               "Lymph_Ecoli_1_25_065",
               "Lymph_Ecoli_1_25_080",
               "Lymph_Ecoli_1_25_082",
               "Lymph_Ecoli_1_25_083",
               "Lymph_Ecoli_1_25_085",
               "Lymph_Ecoli_1_25_087",
               "Lymph_Ecoli_1_25_101",
               "Lymph_Ecoli_1_25_102"
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
