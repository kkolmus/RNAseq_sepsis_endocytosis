### LIBS ###

suppressPackageStartupMessages({
  library(futile.logger)
  library(limma)
  library(Glimma)
  library(edgeR)
  library(Homo.sapiens)
  library(RColorBrewer)
})

### MAIN ###

flog.threshold("DEBUG")

flog.debug("Set directories and read files for analysis")

proj.dir = "~/Desktop/HT projects/RNAseq_sepsis_endocytosis"
raw.data.dir = "~/Desktop/HT projects/RNAseq_sepsis_endocytosis/raw_data"
data.dir = "~/Desktop/HT projects/RNAseq_sepsis_endocytosis/data"

rpm <- readRDS(file.path(data.dir, "RPM.RDS"))
metadata <- readRDS(file.path(data.dir, "metadata.RDS"))


flog.debug("Organising gene annotations")

geneid <- rownames(rpm)
genes <- select(Homo.sapiens, 
                keys = geneid, 
                keytype = "SYMBOL",
                columns = c("ENSEMBL", "ENTREZID", "TXCHROM"))


flog.debug("Organising sample information in metadata")

metadata <- metadata[,c(1,2,48)]

metadata$outcome <- factor(metadata$outcome)
metadata$new_sample_names <- paste0(metadata$geo_accession, "_", metadata$outcome)


flog.debug("Organising sample information in data")

colnames(rpm) <- metadata$new_sample_names 



flog.debug("Prepare separate files")

for (i in names(rpm)) {
  row_names <- rownames(rpm)
  column1 <- rpm[[i]]
  temp_file <- data.frame(GeneIDs = row_names, counts = column1)
  write.table(temp_file, 
              file.path(data.dir, paste0(i, ".txt")), 
              row.names = FALSE,
              quote = FALSE,
              sep = "\t")
}


file.list <- list.files(data.dir, pattern = ".txt")

setwd(data.dir)

# for (f in file.list) {
#   assign(paste0(f), read.table(f))
# }


flog.debug("Read and Merge a Set of Files Containing Count Data")

x <- readDGE(file.list, columns = c(1, 2))
class(x)

colnames(x) <- metadata$new_sample_names
group <- as.factor(metadata$outcome)

x$samples$group <- group


flog.debug("Data pre-processing")

# Transformations from the raw-scale

RPM <- x$counts
class(RPM)

RPM <- as.data.frame(RPM)
class(RPM)

RPM1 <- mutate_all(RPM, funs(.+1))
lRPM <- mutate_all(RPM1, funs(log2))

L <- mean(x$samples$lib.size) * 1e-6
M <- median(x$samples$lib.size) * 1e-6
c(L, M)

summary(lRPM)


flog.debug("Removing genes that are lowly expressed")

table(rowSums(x$counts==0)==129)

keep.exprs <- filterByExpr(x, group = group)
x <- x[keep.exprs, , keep.lib.sizes = FALSE]
dim(x)


