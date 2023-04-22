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

library(iq)
library(stringr)
library(matrixStats)

pathes <- list.files("G:\\Dropbox\\papers_Fengchao\\msfragger_dia\\script\\results\\", pattern = "(report\\.tsv)|(diann-output\\.tsv)", full.names = TRUE, recursive = TRUE)

globalPrecursorFDR <- 0.01
globalProteinFDR <- 0.01
runSpecificPrecursorFDR <- 0.01
runSpecificProteinFDR <- 0.01

for (path in pathes) {
  if (file.exists(path)) {
    if (str_detect(path, "runtime")) {
      next
    }

    out_path <- str_c(dirname(path), "/precursor_maxlfq.tsv")
    if (file.exists(out_path) == FALSE) {
      print(out_path)
      print(str_c("global precursor FDR = ", globalPrecursorFDR))
      print(str_c("global protein FDR = ", globalProteinFDR))
      print(str_c("run specific precursor FDR = ", runSpecificPrecursorFDR))
      print(str_c("run specific protein FDR = ", runSpecificProteinFDR))

      df <- fast_read(path,
                      sample_id = "Run",
                      primary_id = "Precursor.Id",
                      secondary_id = "Precursor.Id",
                      intensity_col = "Precursor.Normalised",
                      annotation_col = NULL,
                      filter_string_equal = NULL,
                      filter_string_not_equal = NULL,
                      filter_double_less = c("Global.Q.Value" = globalPrecursorFDR, "Global.PG.Q.Value" = globalProteinFDR, "Q.Value" = runSpecificPrecursorFDR, "PG.Q.Value" = runSpecificProteinFDR),
                      filter_double_greater = NULL,
                      intensity_col_sep = NULL,
                      intensity_col_id = NULL,
                      na_string = "0")
      df_norm <- fast_preprocess(df$quant_table, median_normalization = FALSE, log2_intensity_cutoff = 0, pdf_out = NULL) # Do not enable median normalization because the tools have their own normalization algorithm.
      df_maxlfq <- fast_MaxLFQ(df_norm, row_names = df$protein[, 1], col_names = df$sample) # MaxLFQ requires a precursor having non-zero intensities in at least two runs, which makes the "quantified" entries slightly fewer than those from the DIA-NN main report.
      maxlfq <- df_maxlfq$estimate
      maxlfq[maxlfq <= 0] <- NA
      write.table(2^maxlfq, out_path, quote = FALSE, sep = "\t", col.names = NA)
    }

    out_path <- str_c(dirname(path), "/modified_sequence_maxlfq.tsv")
    if (file.exists(out_path) == FALSE) {
      print(out_path)
      print(str_c("global precursor FDR = ", globalPrecursorFDR))
      print(str_c("global protein FDR = ", globalProteinFDR))
      print(str_c("run specific precursor FDR = ", runSpecificPrecursorFDR))
      print(str_c("run specific protein FDR = ", runSpecificProteinFDR))

      df <- fast_read(path,
                      sample_id = "Run",
                      primary_id = "Modified.Sequence",
                      secondary_id = "Precursor.Id",
                      intensity_col = "Precursor.Normalised",
                      annotation_col = NULL,
                      filter_string_equal = NULL,
                      filter_string_not_equal = NULL,
                      filter_double_less = c("Global.Q.Value" = globalPrecursorFDR, "Global.PG.Q.Value" = globalProteinFDR, "Q.Value" = runSpecificPrecursorFDR, "PG.Q.Value" = runSpecificProteinFDR),
                      filter_double_greater = NULL,
                      intensity_col_sep = NULL,
                      intensity_col_id = NULL,
                      na_string = "0")
      df_norm <- fast_preprocess(df$quant_table, median_normalization = FALSE, log2_intensity_cutoff = 0, pdf_out = NULL) # Do not enable median normalization because the tools have their own normalization algorithm.
      df_maxlfq <- fast_MaxLFQ(df_norm, row_names = df$protein[, 1], col_names = df$sample) # MaxLFQ requires a precursor having non-zero intensities in at least two runs, which makes the "quantified" entries slightly fewer than those from the DIA-NN main report.
      maxlfq <- df_maxlfq$estimate
      maxlfq[maxlfq <= 0] <- NA
      write.table(2^maxlfq, out_path, quote = FALSE, sep = "\t", col.names = NA)
    }

    out_path <- str_c(dirname(path), "/protein_maxlfq.tsv")
    if (file.exists(out_path) == FALSE) {
      print(out_path)
      print(str_c("global precursor FDR = ", globalPrecursorFDR))
      print(str_c("global protein FDR = ", globalProteinFDR))
      print(str_c("run specific precursor FDR = ", runSpecificPrecursorFDR))
      print(str_c("run specific protein FDR = ", runSpecificProteinFDR))

      df <- fast_read(path,
                      sample_id = "Run",
                      primary_id = "Protein.Group",
                      secondary_id = "Precursor.Id",
                      intensity_col = "Precursor.Normalised",
                      annotation_col = NULL,
                      filter_string_equal = NULL,
                      filter_string_not_equal = NULL,
                      filter_double_less = c("Global.Q.Value" = globalPrecursorFDR, "Global.PG.Q.Value" = globalProteinFDR, "Q.Value" = runSpecificPrecursorFDR, "PG.Q.Value" = runSpecificProteinFDR),
                      filter_double_greater = NULL,
                      intensity_col_sep = NULL,
                      intensity_col_id = NULL,
                      na_string = "0")
      df_norm <- fast_preprocess(df$quant_table, median_normalization = FALSE, pdf_out = NULL) # Do not enable median normalization because the tools have their own normalization algorithm.
      df_maxlfq <- fast_MaxLFQ(df_norm, row_names = df$protein[, 1], col_names = df$sample) # MaxLFQ requires a precursor having non-zero intensities in at least two runs, which makes the "quantified" entries slightly fewer than those from the DIA-NN main report.
      maxlfq <- df_maxlfq$estimate
      maxlfq[maxlfq <= 0] <- NA
      write.table(2^maxlfq, out_path, quote = FALSE, sep = "\t", col.names = NA)
    }

    out_path <- str_c(dirname(path), "/gene_maxlfq.tsv")
    if (file.exists(out_path) == FALSE) {
      print(out_path)
      print(str_c("global precursor FDR = ", globalPrecursorFDR))
      print(str_c("global protein FDR = ", globalProteinFDR))
      print(str_c("run specific precursor FDR = ", runSpecificPrecursorFDR))
      print(str_c("run specific protein FDR = ", runSpecificProteinFDR))

      df <- fast_read(path,
                      sample_id = "Run",
                      primary_id = "Genes",
                      secondary_id = "Precursor.Id",
                      intensity_col = "Precursor.Normalised",
                      annotation_col = NULL,
                      filter_string_equal = NULL,
                      filter_string_not_equal = NULL,
                      filter_double_less = c("Global.Q.Value" = globalPrecursorFDR, "Global.PG.Q.Value" = globalProteinFDR, "Q.Value" = runSpecificPrecursorFDR, "PG.Q.Value" = runSpecificProteinFDR),
                      filter_double_greater = NULL,
                      intensity_col_sep = NULL,
                      intensity_col_id = NULL,
                      na_string = "0")
      df_norm <- fast_preprocess(df$quant_table, median_normalization = FALSE, pdf_out = NULL) # Do not enable median normalization because the tools have their own normalization algorithm.
      df_maxlfq <- fast_MaxLFQ(df_norm, row_names = df$protein[, 1], col_names = df$sample) # MaxLFQ requires a precursor having non-zero intensities in at least two runs, which makes the "quantified" entries slightly fewer than those from the DIA-NN main report.
      maxlfq <- df_maxlfq$estimate
      maxlfq[maxlfq <= 0] <- NA
      write.table(2^maxlfq, out_path, quote = FALSE, sep = "\t", col.names = NA)
    }
  }
}
