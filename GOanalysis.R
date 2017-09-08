library("readr")
library("dplyr")
library("reshape")
library("ggplot2")

# phosphoProtDEforGO <- read_delim("~/Documents/phosphoProt/GOanalysis/phosphoProtDEforGO.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

# read in tab-delimited UniProt results and remove unnecessary columns and column names

uniProtResults <- read_delim("~/Documents/phosphoProt/GOanalysis/uniProtResults.tsv", 
                                  "\t", escape_double = FALSE, trim_ws = TRUE)

# colnames(phosphoProtDEforGO)[1] <- "geneIDs"

# merge the annotation and the differential expression results with left_join
# mega <- left_join(phosphoProtDEforGO, uniProtResults, by = "geneIDs")
# nrow(mega)
# [1] 211

# create GO list for global GO graph and graph the top 70 terms (arbitrarily chosen to make a good looking graph)

GOmega <- unlist(strsplit(mega$`Gene ontology (GO)`, split="; "))

# sort(table(GOmega), decreasing = T)

GOmega <- gsub(" [GO:[0-9]*]", "", GOmega)
grad <- colorRampPalette(c("red", "white"))
par(mar=c(28,4,4,2), cex=.9)
barplot(head(sort(table(GOmega), decreasing = T), n=70), ylim = c(0,50), ylab = "Freqency of GO term", col = grad(70), main = "Top GO terms in Experiment", las=2)

# I wrote out the cleaner version of the uniProt results and metadata from above for later use if needed
# write.table(uniProtResults, file="uniProtResults.tsv", sep='\t', quote=F, row.names = F) 

# write.table(mega, file="GOmetadata.txt", sep='\t',row.names = F, quote=F)

# read in tables listing the GO terms of the genes not present in either situation. Refer to proteomicsAnalysis_part1.R for creating those files.

notInDepleted <- read_csv("~/Documents/phosphoProt/notInDepleted.txt")
notInReplete <- read_csv("~/Documents/phosphoProt/notInReplete.txt")
colnames(notInDepleted)[1] <- "geneIDs"
colnames(notInReplete)[1] <- "geneIDs"

# merge datasheets to produce unique GO terms not found in the opposing samples.

notInDepletedGOterms <- inner_join(notInDepleted, uniProtResults, by="geneIDs")
notInRepleteGOterms <- inner_join(notInReplete, uniProtResults, by="geneIDs")
notInRepleteGOterms <- unique(notInRepleteGOterms)[,] # 24 obs
notInDepletedGOterms <- unique(notInDepletedGOterms)[,] # 14 obs

# subtract the genes present in depleted and replete samples from the Uniprot results, because that listing is of the unique proteins and I have removed the duplicates from the list:

# return all rows from x where there are not matching values in y, keeping just columns from x. Total is 113, as in the uniprot results

shared <- anti_join(uniProtResults, notInDepletedGOterms, by="geneIDs") # 99 obs
shared2 <- anti_join(shared, notInRepleteGOterms, by="geneIDs") # 75 obs
GOForReplete <- rbind(notInDepletedGOterms, shared2)
GOForDepleted <- rbind(notInRepleteGOterms, shared2)

# separate GO terms for each case and sort by most common, starting with DEPLETED

GOmfDsep <- unlist(strsplit(GOForDepleted$`Gene ontology (molecular function)`, split="; "))
GOccDsep <- unlist(strsplit(GOForDepleted$`Gene ontology (cellular component)`, split = "; "))
GObpDsep <- unlist(strsplit(GOForDepleted$`Gene ontology (biological process)`, split = "; "))

topGOccD <- head(sort(table(GOccDsep),decreasing = T), n=8)
topGOmfD <- head(sort(table(GOmfDsep),decreasing = T), n=8)
topGObpD <- head(sort(table(GObpDsep),decreasing = T), n=8)

# visual check to see numbers

topGObpD
# GObpDsep
# positive regulation of transcription elongation from RNA polymerase II promoter [GO:0032968] 3 
# protein phosphorylation [GO:0006468] 3 
# transcription elongation from RNA polymerase II promoter [GO:0006368] 3 
# intracellular signal transduction [GO:0035556] 3
# microtubule-based movement [GO:0007018] 3 
topGOccD
# GOccDsep
# integral component of membrane [GO:0016021]                        nucleus [GO:0005634]                      cytoplasm [GO:0005737] 
# 12                                          10                                           9 
# nucleus [GO:0005634]                plasma membrane [GO:0005886] 
# 8                                           5 
topGOmfD
# GOmfDsep
# ATP binding [GO:0005524]        metal ion binding [GO:0046872]              DNA binding [GO:0003677] 
# 23                                     7                                     6 
# protein kinase activity [GO:0004672]    magnesium ion binding [GO:0000287] 
# 5                                     4 


# same process as above, this time for Replete (control) samples. visual check in there too.
GOmfRsep <- unlist(strsplit(GOForReplete$`Gene ontology (molecular function)`, split="; "))
GOccRsep <- unlist(strsplit(GOForReplete$`Gene ontology (cellular component)`, split="; "))
GObpRsep <- unlist(strsplit(GOForReplete$`Gene ontology (biological process)`, split="; "))

topGOmfR <- head(sort(table(GOmfRsep),decreasing = T), n=8)
topGOccR <- head(sort(table(GOccRsep),decreasing = T), n=8)
topGObpR <- head(sort(table(GObpRsep),decreasing = T), n=8)

topGOmfR
# GOmfRsep
# ATP binding [GO:0005524]            DNA binding [GO:0003677]        kinase activity [GO:0016301] 
# 18                                   5                                   4 
# magnesium ion binding [GO:0000287]      metal ion binding [GO:0046872] 
# 3                                   3 
topGOccR
# GOccRsep
# nucleus [GO:0005634]                      cytoplasm [GO:0005737] integral component of membrane [GO:0016021] 
# 8                                           8                                           8 
# nucleus [GO:0005634]                plasma membrane [GO:0005886] 
# 8                                           5 
topGObpR
# GObpRsep
# positive regulation of transcription elongation from RNA polymerase II promoter [GO:0032968] 3 
# transcription elongation from RNA polymerase II promoter [GO:0006368] 3 
# photosynthesis [GO:0015979] 3 
# DNA replication [GO:0006260] 2 
# embryo development ending in seed dormancy [GO:0009793] 2 

# graph the top GO terms for each of the 3 categories. Object above will not work as is with ggplot so must be melted to reformat the dataframe. also had to add a vector manually with the numbers of gene counts for the depleted and reformatted the GO terms for the x axis labels because some were too long and skewed the margins

bp <- as.data.frame(topGObpR)
bp$Depleted <- c(0, 3, 3, 3, 2, 2, 0, 0)
colnames(bp) <- c("GO term", "Replete", "Depleted")
bp <- melt(bp, id="GO term")
bp$`GO term` <- c("photosynthesis [GO:0015979]", "pos reg of trans elongation from RNA polII prom [GO:0032968]","TOR signaling [GO:0031929]",   "transcription elongation from RNA polII promoter [GO:0006368]","DNA replication [GO:0006260]","electron transport chain [GO:0022900]","embryo development ending in seed dormancy [GO:0009793]","endocytosis [GO:0006897]","photosynthesis [GO:0015979]","pos reg of trans elongation from RNA polII prom [GO:0032968]","TOR signaling [GO:0031929]","transcription elongation from RNA polII promoter [GO:0006368]","DNA replication [GO:0006260]", "electron transport chain [GO:0022900]","embryo development ending in seed dormancy [GO:0009793]","endocytosis [GO:0006897]")
go <- factor(bp$`GO term`)
Condition <- bp$variable
coun <- bp$value
BP <- ggplot(bp, aes(go, coun)) + geom_bar(aes(fill = Condition), width = 0.5, position = position_dodge(width=0.5), stat="identity") + xlab("GO term") + ggtitle("Biological Process") + ylab("Counts") + theme(axis.text.x = element_text(angle=70, color="black", face = "bold", hjust=1), legend.title = element_text(face="bold", size=12), plot.title=element_text(hjust=0.5, size=16, face="bold"), axis.title.x = element_text(face="bold", size=12), axis.title.y= element_text(face="bold", size=12))

# same as above for molecular function, but did not have to rename x axis variables

mf <- as.data.frame(topGOmfR)
mf$Deplete <- c(26, 7, 5, 8, 4, 0, 5, 0)
colnames(mf)<- c("GO term", "Replete", "Depleted")
mf <- melt(mf, id="GO term")
go <- factor(mf$`GO term`)
Condition <- mf$variable
coun <- mf$value
MF <- ggplot(mf, aes(go, coun)) + geom_bar(aes(fill = Condition), width = 0.5, position = position_dodge(width=0.5), stat="identity") + xlab("GO term") + ggtitle("Molecular Function") + ylab("Counts") + theme(axis.text.x = element_text(angle=60, color="black", face = "bold", hjust=1), legend.title = element_text(face="bold", size=12), plot.title=element_text(hjust=0.5, size=16, face="bold"), axis.title.x = element_text(face="bold", size=12), axis.title.y= element_text(face="bold", size=12))

# same as above for cellular component GO terms

cc <- as.data.frame(topGOccR)
cc$Deplete <- c(18,16,11,5,5,5,4,0)
colnames(cc) <- c("GO term", "Replete", "Depleted")
cc <- melt(cc, id="GO term")
go <- factor(cc$`GO term`)
Condition <- cc$variable
coun <- cc$value
CC <- ggplot(cc, aes(go, coun)) + geom_bar(aes(fill = Condition), width = 0.5, position = position_dodge(width=0.5), stat="identity") + xlab("GO term") + ggtitle("Cellular Component") + ylab("Counts") + theme(axis.text.x = element_text(angle=60, color="black", face = "bold", hjust=1), legend.title = element_text(face="bold", size=12), plot.title=element_text(hjust=0.5, size=16, face="bold"), axis.title.x = element_text(face="bold", size=12), axis.title.y= element_text(face="bold", size=12))

##############################################
# Differential Expression GO tables
##############################################

# write out tables with the top DE genes up and down regulated genes in direction control relative to experimental 

uniProtResults <- read_delim("~/Documents/phosphoProt/GOanalysis/uniProtResults.tsv","\t", escape_double = FALSE, trim_ws = TRUE)
phosphoProtDEforGO <- read_delim("~/Documents/phosphoProt/GOanalysis/phosphoProtDEforGO.txt","\t", escape_double = FALSE, trim_ws = TRUE)
p05 <- head(phosphoProtDEforGO, n=54)
p05$logFC <- as.numeric(as.character(p05$logFC))
positive <- p05[p05[2] > 0,] # select positive values = down in experimental
negative <- p05[p05[2] < 0,] # select negative values = up in experimental
colnames(positive)[1] <- "geneIDs"
colnames(negative)[1] <- "geneIDs"
upReg <- left_join(positive, uniProtResults, by="geneIDs")
downReg <- left_join(negative, uniProtResults, by="geneIDs")
write.table(upReg, file="upRegulatedInControl.txt", sep='\t', quote=F)
write.table(downReg, file="downRegulatedInControl.txt", sep='\t', quote=F)