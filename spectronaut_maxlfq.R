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
library(matrixStats)
library(tidyverse)


generate_temp_file <- function(input_path, output_path) {
  xx <- read_tsv(input_path, na = c("NaN", "NA", ""), col_select = c(R.FileName, PG.ProteinGroups, EG.Qvalue, EG.ModifiedSequence, FG.Charge, PG.Qvalue, "PG.QValue (Run-Wise)", "EG.TotalQuantity (Settings)"))
  xx <- xx %>%
    mutate(EG.PrecursorId = paste0(EG.ModifiedSequence, FG.Charge))
  write_tsv(xx, output_path, na = "NA", quote = "none", escape = "none")
}


pathes <- c("G:\\Dropbox\\papers_Fengchao\\msfragger_dia\\script\\results\\benchmark\\spectronaut\\17\\20230227_180032_SN17_RealDilutionSeries_Report.tsv",
            "G:\\Dropbox\\papers_Fengchao\\msfragger_dia\\script\\results\\benchmark\\spectronaut\\14\\20230227_181426_directDIA_LymphEcoli_Report.tsv",
            "G:\\Dropbox\\papers_Fengchao\\msfragger_dia\\script\\results\\2018-HeLa\\spectronaut\\20230315_131200_amodei_triplicate_Report.tsv",
            "G:\\Dropbox\\papers_Fengchao\\msfragger_dia\\script\\results\\2020-Yeast\\spectronaut\\20230315_131543_Searle_2020_quadruplicate_Report.tsv"
)

globalPrecursorFDR <- 0.01
globalProteinFDR <- 0.01
runSpecificProteinFDR <- 0.01

for (path in pathes) {
  if (file.exists(path)) {
    out_path <- str_c(dirname(path), "/precursor_maxlfq.tsv")
    print(out_path)
    print(str_c("global precursor FDR = ", globalPrecursorFDR))
    print(str_c("global protein FDR = ", globalProteinFDR))
    print(str_c("run specific protein FDR = ", runSpecificProteinFDR))

    temp_path <- str_c(path, "_temp")
    generate_temp_file(path, temp_path)

    df <- fast_read(temp_path,
                    sample_id = "R.FileName",
                    primary_id = "EG.PrecursorId",
                    secondary_id = "EG.PrecursorId",
                    intensity_col = "EG.TotalQuantity (Settings)",
                    annotation_col = NULL,
                    filter_string_equal = NULL,
                    filter_string_not_equal = NULL,
                    filter_double_less = c("PG.Qvalue" = globalProteinFDR, "PG.QValue (Run-Wise)" = runSpecificProteinFDR, "EG.Qvalue" = globalPrecursorFDR),
                    filter_double_greater = c("EG.TotalQuantity (Settings)" = 2), # There are weird intensities with value 1, which makes no sense.
                    intensity_col_sep = NULL,
                    intensity_col_id = NULL,
                    na_string = "NaN")
    df_norm <- fast_preprocess(df$quant_table, median_normalization = FALSE, log2_intensity_cutoff = 0, pdf_out = NULL) # Do not enable median normalization because the tools have their own normalization algorithm.
    df_maxlfq <- fast_MaxLFQ(df_norm, row_names = df$protein[, 1], col_names = df$sample) # MaxLFQ requires a precursor having non-zero intensities in at least two runs, which makes the "quantified" entries slightly fewer than those from the Spectronaut main report.
    maxlfq <- df_maxlfq$estimate
    maxlfq[maxlfq <= 0] <- NA
    write.table(2^maxlfq, out_path, quote = FALSE, sep = "\t", col.names = NA)

    out_path <- str_c(dirname(path), "/modified_sequence_maxlfq.tsv")
    print(out_path)
    print(str_c("global precursor FDR = ", globalPrecursorFDR))
    print(str_c("global protein FDR = ", globalProteinFDR))
    print(str_c("run specific protein FDR = ", runSpecificProteinFDR))

    temp_path <- str_c(path, "_temp")
    generate_temp_file(path, temp_path)

    df <- fast_read(temp_path,
                    sample_id = "R.FileName",
                    primary_id = "EG.ModifiedSequence",
                    secondary_id = "EG.PrecursorId",
                    intensity_col = "EG.TotalQuantity (Settings)",
                    annotation_col = NULL,
                    filter_string_equal = NULL,
                    filter_string_not_equal = NULL,
                    filter_double_less = c("PG.Qvalue" = globalProteinFDR, "PG.QValue (Run-Wise)" = runSpecificProteinFDR, "EG.Qvalue" = globalPrecursorFDR),
                    filter_double_greater = c("EG.TotalQuantity (Settings)" = 2), # There are weird intensities with value 1, which makes no sense.
                    intensity_col_sep = NULL,
                    intensity_col_id = NULL,
                    na_string = "NaN")
    df_norm <- fast_preprocess(df$quant_table, median_normalization = FALSE, log2_intensity_cutoff = 0, pdf_out = NULL) # Do not enable median normalization because the tools have their own normalization algorithm.
    df_maxlfq <- fast_MaxLFQ(df_norm, row_names = df$protein[, 1], col_names = df$sample) # MaxLFQ requires a precursor having non-zero intensities in at least two runs, which makes the "quantified" entries slightly fewer than those from the Spectronaut main report.
    maxlfq <- df_maxlfq$estimate
    maxlfq[maxlfq <= 0] <- NA
    write.table(2^maxlfq, out_path, quote = FALSE, sep = "\t", col.names = NA)

    out_path <- str_c(dirname(path), "/protein_maxlfq.tsv")
    print(out_path)
    print(str_c("global precursor FDR = ", globalPrecursorFDR))
    print(str_c("global protein FDR = ", globalProteinFDR))
    print(str_c("run specific protein FDR = ", runSpecificProteinFDR))

    temp_path <- str_c(path, "_temp")
    generate_temp_file(path, temp_path)

    df <- fast_read(temp_path,
                    sample_id = "R.FileName",
                    primary_id = "PG.ProteinGroups",
                    secondary_id = "EG.PrecursorId",
                    intensity_col = "EG.TotalQuantity (Settings)",
                    annotation_col = NULL,
                    filter_string_equal = NULL,
                    filter_string_not_equal = NULL,
                    filter_double_less = c("PG.Qvalue" = globalProteinFDR, "PG.QValue (Run-Wise)" = runSpecificProteinFDR, "EG.Qvalue" = globalPrecursorFDR),
                    filter_double_greater = c("EG.TotalQuantity (Settings)" = 2), # There are weird intensities with value 1, which makes no sense.
                    intensity_col_sep = NULL,
                    intensity_col_id = NULL,
                    na_string = "NaN")
    df_norm <- fast_preprocess(df$quant_table, median_normalization = FALSE, pdf_out = NULL) # Do not enable median normalization because the tools have their own normalization algorithm.
    df_maxlfq <- fast_MaxLFQ(df_norm, row_names = df$protein[, 1], col_names = df$sample) # MaxLFQ requires a precursor having non-zero intensities in at least two runs, which makes the "quantified" entries slightly fewer than those from the Spectronaut main report.
    maxlfq <- df_maxlfq$estimate
    maxlfq[maxlfq <= 0] <- NA
    write.table(2^maxlfq, out_path, quote = FALSE, sep = "\t", col.names = NA)
  }
}
