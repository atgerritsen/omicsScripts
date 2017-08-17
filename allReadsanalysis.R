library("edgeR")
library("BiocGenerics", lib.loc="/opt/modules/devel/R/3.3.2/lib64/R/library")
library("Biobase", lib.loc="/opt/modules/devel/R/3.3.2/lib64/R/library")

#counts <- read.delim("allReadAssembly_counts.tsv", header=T)
counts <- read.delim("ferret_counts.tsv", header=T)
count.mat <- counts[,2:10]
names <- counts$Transcript
y <- DGEList(counts=count.mat, group=c("diapause","diapause", "diapause","estrus", "estrus", "estrus", "pregnant","pregnant","pregnant"), genes=names)

keep <- rowSums(cpm(y)>1) >= 2
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
y$samples
#             group lib.size norm.factors
#Diapause1 diapause 27645996    1.0348445
#Diapause2 diapause 27735506    1.1145644
#Diapause3 diapause 30719415    1.0346491
#Estrus1     estrus 29259443    1.0064111
#Estrus2     estrus 32996514    0.9973108
#Estrus3     estrus 33292257    1.0062075
#Pregnant1 pregnant 25912161    0.8281072
#Pregnant2 pregnant 27260785    1.0829385
#Pregnant3 pregnant 29791449    0.9252152


# make MDS to see how all samples cluster
plotMDS(y, main = "Samples to AllReads")

design <- model.matrix(~0+group, data=y$samples)
colnames(design) <- levels(y$samples$group)
design
# diapause estrus pregnant
# Diapause1        1      0        0
# Diapause2        1      0        0
# Diapause3        1      0        0
# Estrus1          0      1        0
# Estrus2          0      1        0
# Estrus3          0      1        0
# Pregnant1        0      0        1
# Pregnant2        0      0        1
# Pregnant3        0      0        1

levels(y$samples$group)
#[1] "diapause" "estrus"   "pregnant"

# fit to a GLM, because there is more than 1 factor in this experiment
y <- estimateDisp(y, design)
fit <- glmFit(y, design)

# design contrasts
my.contrasts <- makeContrasts(EstrusVDiapause=estrus-diapause, EstrusVPregnant=estrus-pregnant, DiapauseVPregnant=diapause-pregnant, levels=design)
lrt.EstrusVDiapause <- glmLRT(fit, contrast=my.contrasts[,"EstrusVDiapause"])
topTags(lrt.EstrusVDiapause)
# Coefficient:  -1*diapause 1*estrus 
# genes     logFC     logCPM       LR       PValue          FDR
# 1221                CL11198Contig2_1  4.572216  2.9727779 78.38919 8.461378e-19 3.381590e-14
# 59749 k55_Locus_47272_Transcript_1_1 -5.099257  1.0841945 74.26831 6.818951e-18 1.362597e-13
# 6086               CL168916Contig1_1 -6.980790  5.8529275 66.49636 3.505412e-16 4.669793e-12
# 72939 k69_Locus_99038_Transcript_1_1  3.966098  5.7059465 65.42948 6.023126e-16 6.017856e-12
# 49913 k43_Locus_9747_Transcript_61_2  2.719133  4.2940933 64.60580 9.148574e-16 7.312455e-12
# 9661                 CL1Contig1303_1  3.296437  4.3700339 61.87884 3.652532e-15 2.432891e-11
# 60784 k55_Locus_84230_Transcript_8_2  9.482563  2.4185683 60.30521 8.123199e-15 4.637766e-11
# 44250  k37_Locus_808_Transcript_32_3 -6.874959  7.2066023 60.01629 9.407574e-15 4.699671e-11
# 64238  k61_Locus_5770_Transcript_1_1 -4.023603  6.0959528 58.40028 2.138607e-14 9.496602e-11
# 41725 k37_Locus_2837_Transcript_14_3  7.898609 -0.5922017 56.19595 6.559644e-14 2.621562e-10

lrt.EstrusVPregnant <- glmLRT(fit, contrast=my.contrasts[,"EstrusVPregnant"])
lrt.DiapauseVPregnant <- glmLRT(fit, contrast=my.contrasts[,"DiapauseVPregnant"])

evd <- topTags(lrt.EstrusVDiapause, n=nrow(y))
evp <- topTags(lrt.EstrusVPregnant, n=nrow(y))
dvp <- topTags(lrt.DiapauseVPregnant, n=nrow(y))

head(evd$table)
# genes     logFC   logCPM       LR       PValue          FDR
# 1221                CL11198Contig2_1  4.572216 2.972778 78.38919 8.461378e-19 3.381590e-14
# 59749 k55_Locus_47272_Transcript_1_1 -5.099257 1.084194 74.26831 6.818951e-18 1.362597e-13
# 6086               CL168916Contig1_1 -6.980790 5.852927 66.49636 3.505412e-16 4.669793e-12
# 72939 k69_Locus_99038_Transcript_1_1  3.966098 5.705946 65.42948 6.023126e-16 6.017856e-12
# 49913 k43_Locus_9747_Transcript_61_2  2.719133 4.294093 64.60580 9.148574e-16 7.312455e-12
# 9661                 CL1Contig1303_1  3.296437 4.370034 61.87884 3.652532e-15 2.432891e-11

# separate out by false discovery rate 
lrt.evd <- decideTestsDGE(lrt.EstrusVDiapause, adjust.method = "BH", p.value= 0.05)
summary(lrt.evd)
#     [,1] 
# -1   381
# 0  39309
# 1    275

#get rownames and DE genes, plot as MA plot
lrt.evttags <- rownames(y)[as.logical(lrt.evd)]
plotSmear(lrt.EstrusVDiapause, de.tags = lrt.evttags, main = "Estrus Minus Diapause MA")

lrt.evp <- decideTestsDGE(lrt.EstrusVPregnant, adjust.method = "BH", p.value= 0.05)
summary(lrt.evp)
#     [,1] 
# -1  1254
# 0  38035
# 1    676

lrt.evptags <- rownames(y)[as.logical(lrt.evp)]
plotSmear(lrt.EstrusVPregnant, de.tags = lrt.evptags, main = "Estrus Minus Pregnancy MA")

lrt.dvp <- decideTestsDGE(lrt.DiapauseVPregnant, adjust.method = "BH", p.value= 0.05)
summary(lrt.dvp)
#     [,1] 
# -1  1099
# 0  38487
# 1    379
lrt.dvptags <- rownames(y)[as.logical(lrt.dvp)]
plotSmear(lrt.DiapauseVPregnant, de.tags = lrt.dvptags, main = "Diapause Minus Pregnancy MA")

# write out all three pairwise comparisons
# write.table(evd$table, file="EstrusMinusDiapause.txt", sep="\t", quote=F, row.names = F)
# write.table(evp$table, file="EstrusMinusPregnancy.txt", sep="\t", quote=F, row.names = F)
# write.table(dvp$table, file="DiapauseMinusPregnancy.txt", sep="\t", quote=F, row.names = F)
