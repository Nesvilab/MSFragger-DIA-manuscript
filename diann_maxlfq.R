rm(list=ls())

library(diann)
library(stringr)
library(matrixStats)

pathes <- list.files("Z:\\yufe\\results\\msfragger_dia_paper\\", pattern = "(report\\.tsv)|(diann-output\\.tsv)", full.names = TRUE, recursive = TRUE)

globalPrecursorFDR <- 0.01
globalProteinFDR <- 0.01
runSpecificPrecursorFDR <- 0.01
runSpecificProteinFDR <- 0.01

median_normalization <- function (df) {
  medians <- colMedians(df, na.rm = TRUE)
  allMedian <- median(medians, na.rm = TRUE)
  return(allMedian * t(t(df) / medians))
}

for (path in pathes) {
  if (file.exists(path)) {
    if (str_detect(path, "benchmark") || str_detect(path, "runtime")) {
      next
    }

    globalPrecursorFDR <- 0.01
    globalProteinFDR <- 0.01
    runSpecificPrecursorFDR <- 0.01
    runSpecificProteinFDR <- 0.01

    if (str_detect(path, "ccRCC") == FALSE) {
      out_path <- str_c(dirname(path), "/modified_sequence_maxlfq.tsv")
      if (file.exists(out_path) == FALSE) {
        print(out_path)
        print(str_c("global precursor FDR = ", globalPrecursorFDR))
        print(str_c("global protein FDR = ", globalProteinFDR))
        print(str_c("run specific precursor FDR = ", runSpecificPrecursorFDR))
        print(str_c("run specific protein FDR = ", runSpecificProteinFDR))

        df <- diann_load(path)
        df_maxlfq <- diann_maxlfq(df[df$Global.Q.Value < globalPrecursorFDR & df$Global.PG.Q.Value < globalProteinFDR & df$Q.Value < runSpecificPrecursorFDR & df$PG.Q.Value < runSpecificProteinFDR, ], sample.header = "Run", group.header = "Modified.Sequence", id.header = "Precursor.Id", quantity.header = "Precursor.Normalised")
        write.table(df_maxlfq, str_c(dirname(path), "/modified_sequence_maxlfq.tsv"), quote = FALSE, sep = "\t", col.names = NA)
      }

      out_path <- str_c(dirname(path), "/protein_maxlfq.tsv")
      if (file.exists(out_path) == FALSE) {
        print(out_path)
        print(str_c("global precursor FDR = ", globalPrecursorFDR))
        print(str_c("global protein FDR = ", globalProteinFDR))
        print(str_c("run specific precursor FDR = ", runSpecificPrecursorFDR))
        print(str_c("run specific protein FDR = ", runSpecificProteinFDR))

        df <- diann_load(path)
        df_maxlfq <- diann_maxlfq(df[df$Global.Q.Value < globalPrecursorFDR & df$Global.PG.Q.Value < globalProteinFDR & df$Q.Value < runSpecificPrecursorFDR & df$PG.Q.Value < runSpecificProteinFDR, ], sample.header = "Run", group.header = "Protein.Group", id.header = "Precursor.Id", quantity.header = "Precursor.Normalised")
        write.table(df_maxlfq, str_c(dirname(path), "/protein_maxlfq.tsv"), quote = FALSE, sep = "\t", col.names = NA)
      }
    } else {
      out_path_1 <- str_c(dirname(path), "/gene_maxlfq.tsv")
      if (file.exists(out_path_1) == FALSE) {
        print(out_path_1)
        print(str_c("global precursor FDR = ", globalPrecursorFDR))
        print(str_c("global protein FDR = ", globalProteinFDR))
        print(str_c("run specific precursor FDR = ", runSpecificPrecursorFDR))

        df <- diann_load(path)
        df_maxlfq <- diann_maxlfq(df[df$Global.Q.Value < globalPrecursorFDR & df$Global.PG.Q.Value < globalProteinFDR & df$Q.Value < runSpecificPrecursorFDR, ], sample.header = "Run", group.header = "Genes", id.header = "Precursor.Id", quantity.header = "Precursor.Normalised")

        write.table(df_maxlfq, out_path_1, quote = FALSE, sep = "\t", col.names = NA)
      }
    }
  }
}
