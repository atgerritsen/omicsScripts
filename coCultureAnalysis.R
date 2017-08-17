# R version 3.4.1 (2017-06-30)

library("edgeR")
library("BiocGenerics")
library("Biobase")

packageVersion("edgeR")
# [1] ‘3.18.1’
packageVersion("BiocGenerics")
# [1] ‘0.22.0’

# merge kallisto outputs by sample type, starting with cyanobacteria
# cyanoCoCulture4.tsv
# cyanoCoCulture5.tsv
# cyanoCoCulture6.tsv
# cyanoControl1.tsv
# cyanoControl2.tsv
# cyanoControl3.tsv

cyanoControl1 <- read_delim("~/Documents/Omics/CoCultureEcoliCyano/analysis/cyanoControl1.tsv",
                            "\t", escape_double = FALSE, trim_ws = TRUE)
cyanoControl2 <- read_delim("~/Documents/Omics/CoCultureEcoliCyano/analysis/cyanoControl2.tsv",
                            "\t", escape_double = FALSE, trim_ws = TRUE)
cyanoControl3 <- read_delim("~/Documents/Omics/CoCultureEcoliCyano/analysis/cyanoControl3.tsv",
                            "\t", escape_double = FALSE, trim_ws = TRUE)
cyanoCoCulture4 <- read_delim("~/Documents/Omics/CoCultureEcoliCyano/analysis/cyanoCoCulture4.tsv",
                            "\t", escape_double = FALSE, trim_ws = TRUE)
cyanoCoCulture5 <- read_delim("~/Documents/Omics/CoCultureEcoliCyano/analysis/cyanoCoCulture5.tsv",
                              "\t", escape_double = FALSE, trim_ws = TRUE)
cyanoCoCulture6 <- read_delim("~/Documents/Omics/CoCultureEcoliCyano/analysis/cyanoCoCulture6.tsv",
                              "\t", escape_double = FALSE, trim_ws = TRUE)

# remove unnecessary columns, keeping only gene names to merge and count column

gene_names <- cyanoControl1$target_id
cyanoControl1counts <- cyanoControl1$est_counts
cyanoControl2counts <- cyanoControl2$est_counts
cyanoControl3counts <- cyanoControl3$est_counts
cyanoCoCulture4counts <- cyanoCoCulture4$est_counts
cyanoCoCulture5counts <- cyanoCoCulture5$est_counts
cyanoCoCulture6counts <- cyanoCoCulture6$est_counts
cyanoCounts <- cbind(gene_names, cyanoControl1counts, cyanoControl2counts, cyanoControl3counts, 
                       cyanoCoCulture4counts, cyanoCoCulture5counts, cyanoCoCulture6counts)

# quick check to make sure all gene names are in the same order (use "merge" if not)
gene_names2 <- cyanoControl2$target_id
table(gene_names == gene_names2)
# TRUE 
# 2613
table(gene_names == gene_names5)
# TRUE 
# 2613

#write out full counts table for future reference
write.table(cyanoCounts, file="cyanoCountTable.tsv", sep='\t', row.names=F, col.names = T, quote=F)

#before reading this file back in removed the "lcl|" from the gene name column with sed via
# sed -i "" 's/^lcl\|//g' cyanoCountTable.tsv

counts <- read.delim("cyanoCountTable.tsv", header=T)
count.mat <- counts[,2:7]
names <- counts$gene_names
groups <- c("control","control","control","coCulture","coCulture","coCulture")
y <- DGEList(counts=count.mat, group=groups, genes=names)

keep <- rowSums(cpm(y)>1) >= 2
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
y$samples
#                           group lib.size norm.factors
# cyanoControl1counts     control 10043240    1.1304427
# cyanoControl2counts     control  9576281    0.9311429
# cyanoControl3counts     control 10326168    1.0047709
# cyanoCoCulture4counts coCulture  8107088    0.9663966
# cyanoCoCulture5counts coCulture  9041170    1.0025318
# cyanoCoCulture6counts coCulture  8980010    0.9759206

# make MDS to see how all samples cluster
plotMDS(y, main = "Samples to AllReads")

design <- model.matrix(~0+group, data=y$samples)
colnames(design) <- levels(y$samples$group)
design
#                         coCulture control
# cyanoControl1counts           0       1
# cyanoControl2counts           0       1
# cyanoControl3counts           0       1
# cyanoCoCulture4counts         1       0
# cyanoCoCulture5counts         1       0
# cyanoCoCulture6counts         1       0

levels(y$samples$group)
# "coCulture" "control"

# only 1 factor in the experiment, so we can use an exact test for a negative binomial distribution
y <- estimateDisp(y, design)
#fit <- glmFit(y, design)
et <- exactTest(y)
topTags(et)
# Comparison of groups:  control-coCulture 
#                                 genes      logFC   logCPM       PValue          FDR
# 2172 CP000100.1_cds_Synpcc7942_2174_2172  1.3897476 7.131686 3.208715e-19 8.377955e-16
# 2170 CP000100.1_cds_Synpcc7942_2172_2170  1.1862685 5.301177 2.840354e-11 3.708083e-08
# 1647 CP000100.1_cds_Synpcc7942_1649_1647 -0.9325372 5.992653 8.444717e-11 7.349719e-08
# 2171 CP000100.1_cds_Synpcc7942_2173_2171  1.1510759 6.751913 5.833474e-09 3.579808e-06
# 1866 CP000100.1_cds_Synpcc7942_1868_1866 -0.7989046 4.346062 6.855243e-09 3.579808e-06
# 1463 CP000100.1_cds_Synpcc7942_1465_1463 -0.6094023 6.951006 2.781684e-07 1.210496e-04
# 19     CP000100.1_cds_Synpcc7942_0019_19 -0.7349507 8.085330 3.715094e-07 1.385730e-04
# 840   CP000100.1_cds_Synpcc7942_0841_840 -0.6139898 6.670876 3.230177e-06 1.007180e-03
# 2353 CP000100.1_cds_Synpcc7942_2355_2353 -0.7448914 7.079053 3.815610e-06 1.007180e-03
# 2570 CP000100.1_cds_Synpcc7942_2572_2570 -1.4938120 6.729360 3.857448e-06 1.007180e-03

cyanoTable <- topTags(et, n=nrow(y))
head(cyanoTable$table)
# output is same as above

# write out full differential expression table
write.table(cyanoTable$table, file="CyanoComparison.txt", sep='\t', quote=F, row.names=F)

# separate out by false discovery rate 
fdr <- decideTestsDGE(et, adjust.method = "BH", p.value= 0.05)
summary(fdr)
#     coCulture+control
# -1                39
# 0               2557
# 1                 15

# get rownames and DE genes, plot as MA plot
fdr.tags <- rownames(y)[as.logical(fdr)]
ex <- plotSmear(et, de.tags = fdr.tags, main = "Cyanobacteria Differential Expression")

# attach annotations post hoc to the diffy exp table: did some sed in terminal to remove leading
# characters and add commas between all brackets to aid in column separation. also need to read 
# original file back in.

library(readr)
selongatusAnnots <- read_csv("~/Documents/Omics/CoCultureEcoliCyano/references/selongatusAnnots.txt", 
                             col_names = FALSE)
colnames(selongatusAnnots) <- c("genes", "ID", "protein", "ABB", "loc", "loc2")
merged <- merge(CyanoComparison, selongatusAnnots, by="genes")
mo <- merged[order(merged$FDR),]
write.table(mo, file="cyanoAnnotatedComparison.txt", sep='\t', quote=F, row.names=F)
