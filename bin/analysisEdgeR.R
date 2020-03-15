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

genes <- genes[!duplicated(genes$SYMBOL),]


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

rm(column1, i, row_names, temp_file)

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

x$samples
dim(x)

flog.debug("Data pre-processing")

# Add gene information

x$genes <- genes

# Transformations from the raw-scale

RPM <- x$counts
class(RPM)

RPM <- as.data.frame(RPM)
class(RPM)

RPM1 <- mutate_all(RPM, funs(.+1))
lRPM <- mutate_all(RPM1, funs(log2))
rownames(lRPM) <- genes$SYMBOL

L <- mean(x$samples$lib.size) * 1e-6
M <- median(x$samples$lib.size) * 1e-6
c(L, M)

summary(lRPM)


flog.debug("Removing genes that are lowly expressed")

table(rowSums(x$counts==0)==129)

keep.exprs <- filterByExpr(x, group = group)
x <- x[keep.exprs, , keep.lib.sizes = FALSE]
dim(x)


flog.debug("Normalizing gene expression distributions")

x <- calcNormFactors(x, method = "TMM")
x$samples$norm.factors
dim(x)

x2 <- x
x2$samples$norm.factors <- 1
x2$counts[,1] <- ceiling(x2$counts[,1]*0.05)
x2$counts[,2] <- x2$counts[,2]*5

par(mfrow=c(1,2))

lcpm <- cpm(x2, log=TRUE)
boxplot(lcpm, las=2, col=col, main="")
title(main="A. Unnormalized data",ylab="Log-cpm")

x2 <- calcNormFactors(x2)  
x2$samples$norm.factors

lcpm <- cpm(x2, log=TRUE)
boxplot(lcpm, las=2, col=col, main="")
title(main="B. Normalized data",ylab="Log-cpm")


flog.debug("Differential expression analysis")

flog.debug("Creating a design matrix and contrasts")

design <- model.matrix(~0+group)
colnames(design) <- gsub("group", "", colnames(design))
design


contrast.matrix <- makeContrasts(
  SurvivedvsDeath = Sepsis_Survival - Sepsis_Death, 
  SurvivedvsSIRS = Sepsis_Survival - SIRS, 
  SIRSvsDeath = Sepsis_Survival - Sepsis_Death, 
  levels = colnames(design))
contrast.matrix


flog.debug("Removing heteroscedascity from count data")

par(mfrow=c(1,2))
v <- voom(x, design, plot = TRUE)
v

vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts = contrast.matrix)
efit <- eBayes(vfit)

plotSA(efit, main="Final model: Mean-variance trend")

# Means (x-axis) and variances (y-axis) of each gene are plotted to show the dependence 
# between the two before voom is applied to the data (left panel) and how the trend is 
# removed after voom precision weights are applied to the data (right panel).


flog.debug("Identify differentially regulated genes")

summary(decideTests(efit))

tfit <- treat(vfit, lfc=1)
dt <- decideTests(tfit)
summary(dt)

de.common <- which(dt[,1]!=0 & dt[,2]!=0)
length(de.common)

head(tfit$genes$SYMBOL[de.common], n = 25)



### SESSION INFO ###

flog.debug("Session Info")
sessionInfo()