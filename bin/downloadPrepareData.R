### LIBS ###

suppressPackageStartupMessages({
  # Formatting/documentation packages
  library(futile.logger)
  library(dplyr)
  library(readxl)
  # Annotation and data import packages
  library(GEOquery)
})


### MAIN ###

flog.threshold("DEBUG")

flog.debug("Set directories")

proj.dir = "~/Desktop/HT projects/RNAseq_sepsis_endocytosis"
raw.data.dir = "~/Desktop/HT projects/RNAseq_sepsis_endocytosis/raw_data"
data.dir = "~/Desktop/HT projects/RNAseq_sepsis_endocytosis/data"
geo = "GSE63042"


flog.debug("Download expression data")

getGEOSuppFiles(geo, makeDirectory = FALSE, baseDir = raw.data.dir)


flog.debug("Prepare expression data")

gz <- list.files(path = raw.data.dir, pattern = ".gz", full.names = TRUE)
lapply(gz, gunzip)


flog.debug("Load file with RNA-seq data")
temp_file <- list.files(path = raw.data.dir, pattern = ".xlsx", full.names = TRUE)
exprs <- read_excel(file.path(temp_file))
exprs <- exprs[, c(2, 7:ncol(exprs))]

# number of identified reads
unfiltered_df <- nrow(exprs)
# remove unidentified reads
exprs <- exprs[exprs$gene_name != "-",]
# number of reads which map to genes
filtered_df <- nrow(exprs)

flog.debug("Percent of mapped genes")
mapping = round((filtered_df/unfiltered_df) * 100, 2)
print(paste0("Number of reads which map to genes: ", mapping, " %"))

# remove duplicated reads
exprs <- exprs[! duplicated(exprs$gene_name),]
dup_rm <- nrow(exprs)

flog.debug("Number of duplicated genes")
dup = round((100 -(dup_rm/filtered_df)*100), 2)
print(paste0("Number of duplicated reads amongst the identified genes: ", dup, " %"))


flog.debug("Prepare final dataframe")
exprs <- as.data.frame(exprs)
rownames(exprs) <- exprs$gene_name
exprs <- exprs[, -1]


saveRDS(exprs, file.path(data.dir, "RPM.RDS"))


flog.debug("Download metadata")

gse <- getGEO(geo, destdir = raw.data.dir, GSEMatrix = TRUE)


flog.debug("Prepare metadata")

metadata <- gse[["GSE63042_series_matrix.txt.gz"]]@phenoData@data
metadata <- mutate(metadata, outcome = with(metadata, 
                                            ifelse(grepl("Septic shock", metadata$`sirs outcomes:ch1`), "Sepsis_Survival",
                                                   ifelse(grepl("severe sepsis", metadata$`sirs outcomes:ch1`), "Sepsis_Survival",
                                                          ifelse(grepl("Uncomplicated sepsis", metadata$`sirs outcomes:ch1`), "Sepsis_Survival",
                                                                 ifelse(grepl("SIRS", metadata$`sirs outcomes:ch1`), "SIRS",
                                                                        "Sepsis_Death"))))))
rownames(metadata) <- metadata$title
metadata <- metadata[match(colnames(exprs), rownames(metadata)), ]


saveRDS(metadata, file = file.path(data.dir, "metadata.RDS"))


flog.debug("Check alignment between expression and metadata")

stopifnot(rownames(metadata) == colnames(exprs))


### SESSION INFO ###
sessionInfo()