
#################################################################################
# Part 2: DE analysis
#################################################################################

# DESeq mentioned as being appropriate for proteomics comparisons, although data have to be in raw counts form. First steps are identical to part 1 to read in the original table and remove all irrelevant rows. Skipping the log+1 normalization to retain the raw read counts for DESeq.

library("readxl")
phosphoprotAnnotated_EK <- read_excel("~/Documents/phosphoProt/phosphoprotAnnotated_EK.xlsx")

# pare down the sheet to names and intensities

protName <- phosphoprotAnnotated_EK$`match to original proteomics`
d1 <- phosphoprotAnnotated_EK$`Intensity D_1`
d2 <- phosphoprotAnnotated_EK$`Intensity D_2`
d3 <- phosphoprotAnnotated_EK$`Intensity D_3`
r1 <- phosphoprotAnnotated_EK$`Intensity R_1`
r2 <- phosphoprotAnnotated_EK$`Intensity R_2`
r3 <- phosphoprotAnnotated_EK$`Intensity R_3`
intensities <- cbind(protName, d1, d2, d3, r1, r2, r3)
intensities <- as.data.frame(intensities)
head(intensities)
# protName                                                                                d1        d2        d3        r1        r2        r3
# 1                                                                                  nope 478420000 475350000  74781000 239260000 132820000  59604000
# 2                             Top Blastx Hit: >ref|NP_001105858.1| LOC732762 [Zea mays] 369120000 336120000 228890000 454470000 427010000 279500000
# 3 Top Blastx Hit: >ref|XP_002511578.1| RNA binding protein, putative [Ricinus communis] 333780000 291340000 219430000 105580000 146460000  99967000
# 4                             Top Blastx Hit: >ref|NP_001105858.1| LOC732762 [Zea mays] 319960000 308260000 228450000 342310000 359690000 217160000
# 5        Top Blastx Hit: >ref|XP_001692808.1| apoferredoxin [Chlamydomonas reinhardtii] 315210000 350060000 262830000 479610000 519120000 338130000
# 6          Top Blastx Hit: >ref|XP_003074063.1| malic enzyme (ISS) [Ostreococcus tauri] 277670000 172920000  67879000 132030000  68875000  52637000

# subset the dataset to remove and keep the "nope" and "No Blastx Hits" IDs in a separate dataset with corresponding row names 

noidea <- intensities$protName == "No Blastx Hits" | intensities$protName == "nope"
noIDs <- intensities[noidea,]
head(noIDs)
# protName          d1        d2        d3        r1        r2       r3
# 1            nope 478420000 475350000  74781000 239260000 132820000 59604000
# 10           nope 226250000         0         0  87975000         0        0
# 25           nope  98129000  97754000 173070000  95445000 166710000 92042000
# 31           nope  85767000  88880000 162550000  85256000 155020000 85223000
# 42 No Blastx Hits  53359000  33067000  19098000  60522000  18332000 14596000
# 43           nope  52116000  45705000  50345000  26747000  32962000 16658000

write.table(noIDs, file="unIdentifiedProteins.txt", sep='\t', quote=F, col.names = TRUE)

# keep all rows with protein IDs and remove all rows where row sums are 0
intensities <- intensities[!noidea,]
intensities$d1 <- as.numeric(levels(intensities$d1))[intensities$d1]
intensities$d2 <- as.numeric(levels(intensities$d2))[intensities$d2]
intensities$d3 <- as.numeric(levels(intensities$d3))[intensities$d3]
intensities$r1 <- as.numeric(levels(intensities$r1))[intensities$r1]
intensities$r2 <- as.numeric(levels(intensities$r2))[intensities$r2]
intensities$r3 <- as.numeric(levels(intensities$r3))[intensities$r3]
intensities$sum <- rowSums(intensities[,2:7])
keepIns <- na.omit(subset(intensities, intensities[,8] != 0))
head(keepIns)
# protName        d1
# Top Blastx Hit: >ref|NP_001105858.1| LOC732762 [Zea mays]                             369120000
# Top Blastx Hit: >ref|XP_002511578.1| RNA binding protein, putative [Ricinus communis] 333780000
# Top Blastx Hit: >ref|NP_001105858.1| LOC732762 [Zea mays]                             319960000
# Top Blastx Hit: >ref|XP_001692808.1| apoferredoxin [Chlamydomonas reinhardtii]        315210000
# Top Blastx Hit: >ref|XP_003074063.1| malic enzyme (ISS) [Ostreococcus tauri]          277670000
# Top Blastx Hit: >ref|XP_001756574.1| transcription elongation factor SPT6 [Physcomitrella patens subsp. patens] 254490000
# d2        d3        r1        r2        r3        sum
# 2 336120000 228890000 454470000 427010000 279500000 2095110000
# 3 291340000 219430000 105580000 146460000  99967000 1196557000
# 4 308260000 228450000 342310000 359690000 217160000 1775830000
# 5 350060000 262830000 479610000 519120000 338130000 2264960000
# 6 172920000  67879000 132030000  68875000  52637000  772011000
# 7 230270000 184130000 134810000 112180000  57322000  973202000

fors <- keepIns[(keepIns$d1 !=0 & keepIns$d2 !=0 & keepIns$d3 !=0) | (keepIns$d1==0 & keepIns$d2!=0 & keepIns$d3!=0) | (keepIns$d1!=0 & keepIns$d2==0 & keepIns$d3!=0) | (keepIns$d1 !=0 & keepIns$d2 !=0 & keepIns$d3==0) | (keepIns$r1==0 & keepIns$r2==0 & keepIns$r3==0) | (keepIns$d1==0 & keepIns$d2==0 & keepIns$d3 == 0),]
rawCounts <- fors[(fors$r1 !=0 & fors$r2 !=0 & fors$r3 !=0) | (fors$r1==0 & fors$r2!=0 & fors$r3!=0) | (fors$r1!=0 & fors$r2==0 & fors$r3!=0) | (fors$r1!=0 & fors$r2!=0 & fors$r3==0) | (fors$d1==0 & fors$d2==0 & fors$d3==0 & fors$r1 !=0 & fors$r2!=0 & fors$r3 ==0)| (fors$d1==0 & fors$d2==0 & fors$d3==0 & fors$r1 !=0 & fors$r2==0 & fors$r3 !=0)| (fors$d1==0 & fors$d2==0 & fors$d3==0 & fors$r1 ==0 & fors$r2!=0 & fors$r3 !=0) | (fors$r1==0 & fors$r2==0 & fors$r3 == 0 & fors$d1 != 0 & fors$d2 !=0 & fors$d3==0)| (fors$r1==0 & fors$r2==0 & fors$r3 == 0 & fors$d1 == 0 & fors$d2 !=0 & fors$d3!=0)| (fors$r1==0 & fors$r2==0 & fors$r3 == 0 & fors$d1 != 0 & fors$d2 ==0 & fors$d3!=0) | (fors$r1==0 & fors$r2==0 & fors$r3 == 0 & fors$d1 != 0 & fors$d2 !=0 & fors$d3!=0) ,]

nonNormal <- unique(nonNormal[,])

write.table(nonNormal, file="nonNormalizedCountsForDE.txt", sep='\t',quote=F, row.names = F)
countTable <- read.delim("nonNormalizedCountsForDE.txt", header=T, sep='\t')

library("edgeR")

count.mat <- countTable[,3:8]
names <- countTable$X1

groups <- c("depleted","depleted","depleted","replete","replete","replete")
y <- DGEList(counts=count.mat, group=groups, genes=names)

keep <- rowSums(cpm(y)>1) >= 2
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
y$samples
# group   lib.size norm.factors
# d1 depleted 7048155600    1.0550214
# d2 depleted 7130516100    1.1086725
# d3 depleted 8417022500    1.1797161
# r1  replete 6999899400    0.9203917
# r2  replete 7967988300    0.8559535
# r3  replete 4863590200    0.9198883
plotMDS(y, main = "Samples to AllReads")
# design <- model.matrix(~0+group, data=y$samples)
# colnames(design) <- levels(y$samples$group)
# design
# depleted replete
# d1        1       0
# d2        1       0
# d3        1       0
# r1        0       1
# r2        0       1
# r3        0       1

y <- estimateCommonDisp(y)
y <- estimateTagwiseDisp(y)
et <- exactTest(y)
topTags(et)
phosphoProtTable <- topTags(et, n=nrow(y))
write.table(phosphoProtTable$table, file="phosphoProtDEcomparison.txt", sep='\t', quote=F, row.names=F)
fdr <- decideTestsDGE(et, adjust.method = "BH", p.value= 0.05)
summary(fdr)
# depleted+replete
# -1               17
# 0               176
# 1                11