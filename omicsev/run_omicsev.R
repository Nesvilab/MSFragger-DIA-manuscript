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

rm(list=ls())

library(data.table)
library(tidyverse)
library(stringr)
library(OmicsEV)

setwd("/storage/yufe/results/msfragger_dia_paper/ccRCC/omicsev/OmicsEV/files/")

print("Generating sample_list.tsv from the col_annot.tsv")
sl <- fread("col_annot.tsv", stringsAsFactors = F, data.table = F)
samplist <- sl %>%
  dplyr::select(caseID, type) %>%
  dplyr::mutate(batch = 1, order = row_number()) %>%
  dplyr::rename(sample = caseID, class = type)
write.table(samplist, "sample_list.tsv", quote = F, row.names = F, sep = "\t")

dir.create(file.path("datasets2"))

print("Processing TMT table")
al_tmt <- fread("datasets/tmt_abundance_gene_MD_Alexey_20220814.tsv", stringsAsFactors = F, data.table = F)

cols <- data.frame(prot_aliquot = colnames(al_tmt)[-c(1:5)], stringsAsFactors = F) %>%
  dplyr::left_join(sl, by = "prot_aliquot") %>%
  dplyr::select(prot_aliquot, caseID) %>%
  na.omit()

al_tmt <- al_tmt[, c("Index", cols$prot_aliquot)]
colnames(al_tmt) <- c("ID", cols$caseID)
al_tmt[,-1] <- 2^(al_tmt[, -1]) # TMT table's intensity has log2 transferred.

# CPT0081990003 (C3N-00646-T) and CPT0086890003 (C3L-01836-N) are mislabeled. Should be switched.
wt <- which(colnames(al_tmt)=="C3N-00646-T")
wn <- which(colnames(al_tmt)=="C3L-01836-N")
colnames(al_tmt)[wt] <- "C3L-01836-N"
colnames(al_tmt)[wn] <- "C3N-00646-T"

write.table(al_tmt, "datasets2/tmt.tsv", quote = F, row.names = F, sep = "\t")


inputList <- c("datasets/diann_gene_maxlfq.tsv",
               "datasets/diaumpire_gene_maxlfq.tsv",
               "datasets/diaumpiredda_gene_maxlfq.tsv",
               "datasets/msfraggerdia_gene_maxlfq.tsv",
               "datasets/msfraggerdiadda_gene_maxlfq.tsv")

outputList <- c("datasets2/diann.tsv",
                "datasets2/diaumpire.tsv",
                "datasets2/diaumpiredda.tsv",
                "datasets2/msfraggerdia.tsv",
                "datasets2/msfraggerdiadda.tsv")

for (i in seq(1, length(inputList))) {
  print(str_c("Processing ", inputList[i]))

  diaInput <- fread(inputList[i], stringsAsFactors = F, data.table = F)

  colnames(diaInput) <- gsub("CPTAC_CCRCC_W_JHU_20190112_LUMOS_", "", colnames(diaInput))
  colnames(diaInput) <- gsub(".mzML", "", colnames(diaInput))
  colnames(diaInput) <- gsub("_T", "-T", colnames(diaInput))
  colnames(diaInput) <- gsub("_NAT", "-N", colnames(diaInput))
  colnames(diaInput)[1] <- "ID"

  # change the N and T of 00360. They are mislabeled.
  wt <- which(colnames(diaInput)=="C3L-00360-T")
  wn <- which(colnames(diaInput)=="C3L-00360-N")
  colnames(diaInput)[wt] <- "C3L-00360-N"
  colnames(diaInput)[wn] <- "C3L-00360-T"

  # There are row with empty gene names. They are very rare but crash OmicsEV. Remove them.
  empty_rows <- diaInput[, 1] == ""
  diaInput <- diaInput[!empty_rows, ]

  write.table(diaInput, outputList[i], quote = F, row.names = F, sep = "\t")
}

run_omics_evaluation(data_dir = "datasets2/",
                     sample_list = "sample_list.tsv",
                     x2 = "rpkm_rna_20220810.tsv",
                     cpu = 56,
                     method_for_fun = "co",
                     data_type = "gene",
                     class_for_cor = "Tumor",
                     class_for_fun = "Tumor",
                     class_for_ml = NULL,
                     use_existing_data=TRUE)

