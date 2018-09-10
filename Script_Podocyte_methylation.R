library(minfi)
library(IlluminaHumanMethylationEPICmanifest)
library(RColorBrewer)
library(limma)
library("IlluminaHumanMethylationEPICanno.ilm10b2.hg19")
library(wateRmelon)
library(methylumi)
library(methyAnalysis)
library(gplots)
library(sva)
library(ggplot2)
library(Heatplus)
library(methylKit)
library(ENmix)
library(Gviz)

# Set working directory where your IDAT files and samplesheet are located
# Replace back slash with forward slash
setwd("/panfs/panasas01/sscm/sa16666/450kanalysis/Podocyte_Methylation/")

# Place the samplesheet and IDAT files in one folder
# Important: There should be no other file present in this forlder except the samplesheet and IDAT files
# Make that folder as the working directory as follows:
#setwd("~/data/IDAT_1000/sub") # Go to the folder and copy the location, Replace all backslash (\) with forward slash (/)
baseDir <- "/panfs/panasas01/sscm/sa16666/450kanalysis/Podocyte_Methylation/" # Call the same location as basedirectory
list.files(baseDir) # See the IDAT files and samplesheet files listed 

target <- read.metharray.sheet(baseDir) # It will automatically read the samplesheet file saved as .csv
baseDir <- "/panfs/panasas01/sscm/sa16666/450kanalysis/Podocyte_Methylation/" # Call the same location as basedirectory

rgset<-read.metharray.exp(baseDir, force=TRUE)
#target$Basename<-target$barcode # If you see problem with loading files in the next step
#rgset<-read.metharray.exp(target=target) # It will automatically read the IDAt files present in the folder
rgset@annotation=c(array="IlluminaHumanMethylationEPIC",annotation="ilm10b2.hg19")

dim(rgset) # Will tell you about the number of probes (CpG sites) and samples in the data
save(rgset, file="rgset_280618.RData") # Save the file you have loaded 
#load("rgset.RData")
#load("detP_rgset.RData")
# Extract all the information required from the rgset data
detP<-detectionP(rgset)
save(detP, file="detP_rgset_280618.RData")
# Vizualising the setection Pvalues
pal <- brewer.pal(8,"Dark2")
tiff("barplot_PVals.tif", res=300, compression = "lzw", height=5, width=10, units="in")
barplot(colMeans(detP), col=pal[factor(target$Group)], las=2,
        cex.names=0.8,ylab="Mean detection p-values") # Here the Sample_Group refers to "Normal", "OPL" and "Cancer" categories
abline(h=0.01,col="red")
legend("topright", legend=levels(factor(target$Sample_Group)), fill=pal,
       bg="white")
dev.off()
# Generating QC report
qcReport(rgset, sampNames=target$Sample_Name, sampGroups=target$Group,
         pdf="qcReport.pdf")


############ Exclude samples  that performed poorly 
##### Samples with >5% missing CpG sites (detection p-value>0.01) will be excluded
keep <- colMeans(detP) < 0.05
rgSet_1 <- rgset[,keep]
dim(rgSet_1)

# remove poor quality samples from target data
# see if you have lost some samples' If yes for further analysis use target_1 instead of targtes
target_1 <- target[keep,]
dim(target_1)

# remove poor quality samples from detection p-value table
detP_1 <- detP[,keep]
dim(detP_1)

##### Normalize the data###

mset<-preprocessFunnorm(rgSet_1, nPCs=2, sex = NULL, bgCorr = TRUE, dyeCorr = TRUE, verbose = TRUE)

##### Not normalized data fro comparison###

mSetRaw <- preprocessRaw(rgSet_1)

# Do not do this step
# visualise what the data looks like before and after normalisation
#par(mfrow=c(1,2))
#densityPlot(rgSet_1, sampGroups=target_1$Sample_Group,main="Raw", legend=FALSE)
#legend("top", legend = levels(factor(target_1$Sample_Group)),
#       text.col=brewer.pal(8,"Dark2"))
#densityPlot(getBeta(mset), sampGroups=target_1$Sample_Group,main="Normalized", legend=FALSE)
#legend("top", legend = levels(factor(target_1$Sample_Group)),
#      text.col=brewer.pal(8,"Dark2"))

# MDS plots to look at largest sources of variation
#par(mfrow=c(1,2))
# Is it Sample Group ?
#plotMDS(getM(mset), top=1000, gene.selection="common",
#       col=pal[factor(target_1$Sample_Group)])
#legend("top", legend=levels(factor(target_1$Sample_Group)), text.col=pal,bg="white", cex=0.7)

# Is it Sentrix_ID, the array used for the hybridization ?
#plotMDS(getM(mset), top=1000, gene.selection="common",
#       col=pal[factor(target_1$Sentrix_ID)])
#legend("top", legend=levels(factor(target_1$Sentrix_ID)), text.col=pal,bg="white", cex=0.7)

# Is it Pool_ID, the bisulphite plate used for conversion ?
#plotMDS(getM(mset), top=1000, gene.selection="common",
#       col=pal[factor(target_1$Pool_ID)])
#legend("top", legend=levels(factor(target_1$Pool_ID)), text.col=pal,bg="white", cex=0.7)




# Filtering faulty probes

# ensure probes are in the same order in the mset and detP objects
detP_1 <- detP_1[match(featureNames(mset),rownames(detP_1)),]

# remove any probes that have failed in one or more samples
keep <- rowSums(detP_1 < 0.01) == ncol(mset)
table(keep)
mSetSqFlt <- mset[keep,]
mSetSqFlt

# Check if the predicted sex is right
mSetSqFlt@colData$predictedSex

# remove Sex Chromsomes
# if your data includes males and females, remove probes on the sex chromosomes

#annEPIC = getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
#head(annEPIC)
#keep <- !(featureNames(mSetSqFlt) %in% annEPIC$Name[annEPIC$chr %in% c("chrX","chrY")])
#table(keep)
#mSetSqFlt <- mSetSqFlt[keep,]

# Remove SNP containing probes
# remove probes with SNPs at CpG site customize it!
# with minor allele frequency > 5% for South Asian Samples at target CpG sites
# and at single base extension sites of Type I probes (n=772)
mSetSqFlt_SNP <- dropLociWithSnps(mSetSqFlt)
mSetSqFlt
mSetSqFlt_SNP
# or customize it
#SNP <- read.csv("SNP_EUR_race_5percent.csv",header=F)$V1
#as.character(SNP)
#mSetSqFlt_SNP<- mSetSqFlt[ ! featureNames(mSetSqFlt) %in% SNP, ]

# Exclude cross-reactive probes cutomize it!

Cross_reactive <- read.csv("EPIC_Cross_Reactive.csv",header=F)$V1
as.character(Cross_reactive)
mSetSqFlt_CR<- mSetSqFlt_SNP[ ! featureNames(mSetSqFlt_SNP) %in% Cross_reactive, ]

save(mSetSqFlt_CR, file="mSetSqFlt_CR_280618.RData")

# # calculate M-values for statistical analysis
# Load methylation data
load("mSetSqFlt_CR_280618.RData")
meth <- getM(mSetSqFlt_CR)
head(meth[,1:5])
meth[!is.finite(meth)] <- 0
# Load phenotype data
target<-read.csv("Sample_sheet.csv", header=T)


# See possible effects of technical biological factors before running the analysis:
setwd("C:/MY documents/Podocyte_methylation")

# 1. ComBat correct for date of experiment
batch = target$Batch
combat_meth = ComBat(dat=meth, batch=batch, mod=NULL, par.prior=TRUE, prior.plots=FALSE)
mod<-model.matrix(~as.factor(Group), data=target)
fit_sva= lmFit(combat_meth,mod)
contrast.matrix<-c(0,1)
fitContrasts_sva = contrasts.fit(fit_sva,contrast.matrix)
eb_sva = eBayes(fitContrasts_sva)
top_sva<-topTable(eb_sva, adjust="BH",number=1000000,p.value=0.01,genelist=eb_sva$genes)
# Write table with FDR <=0.05
write.table(top_sva, file="top_ComBat_Batch_050718.txt", row.names=T, sep="\t")
# This part of analysis was run on Bluecrystal as there were issues with loading Meffil locally
library(meffil)
EPIC_anno<-meffil.get.features("epic")
save(EPIC_anno, file="EPIC_anno_meffil.RData")

# Load EPIC_Anno data to identify genes and their loaction both top_sva and top2
# For all significant CpGs
load("EPIC_anno_meffil.RData")
rownames(EPIC_anno)<-EPIC_anno$name
# Merge top_sva with EPIC_anno
de <- merge(top_sva, EPIC_anno, by="row.names",all.x= T)
write.table(de, file="top_ComBat_Batch_annotated_050718.txt", row.names=T, sep="\t")

#######################
# Volcano plot
#########################
head(fit_sva$coef[,2])
head(eb_sva$p.value)
library(rafalib)
tiff("Volcanoplot_All_090718.tif", res=300, compression = "lzw", height=5, width=10, units="in")
splot(fit_sva$coef[,2],-log10(eb_sva$p.value[,1]),col=rgb(0,0,0,.25),xlab="Effect size",ylab="-log10 p-value")
dev.off()

####################
# QQ plot for the findings
######################

library(qqman)
tiff("qqplot_ComBat_Batch.tif", res=300, compression = "lzw", height=5, width=10, units="in")
qq(eb_sva$p.value)
dev.off()

# To get betas from Meth values
betas<-ilogit2(meth)
rownames(betas)
betas<-as.data.frame(betas)

Sel<-betas[(rownames(betas)%in%rownames(top_sva)),]
tbetas<-data.frame(t(Sel)) # transpose betas and merge with pData 
rownames(target)<-target$Basename
Table <- cbind(target,tbetas[match(rownames(target),rownames(tbetas)),])
means<-by(Table[,9:11803], Table$Group, colMeans)
means<-as.table(means)
# use the threshold that you want. In this case 20% 
cond <-(means$"IR"-means$"IS">=0.2)|(means$"IR"-means$"IS"<=-0.2)
count(cond) # 256 sites
cond

# Select the CpGs which have a delta betas of more than 20%
# how many sites fill this condition? : 256 
cond_new<-names(cond)[cond=="TRUE"] 
# top2 has the DMPs selected for FDR and deltabeta 20%
top2<-top_sva[rownames(top_sva)%in%cond_new,]
write.table(top2, file="top2_ComBat_Batch_deltaBetas20_050718.txt", row.names=T, sep="\t")

# get the data from these 256 sites from Beta values 
top2_Sel<-betas[(rownames(betas)%in%rownames(top2)),]

tbetas_top2_Sel<-data.frame(t(top2_Sel)) # transpose betas and merge with pData 
rownames(target)<-target$Basename
Table_top2_Sel <- cbind(target,tbetas_top2_Sel[match(rownames(target),rownames(tbetas_top2_Sel)),])
save(Table_top2_Sel, file="Table_top2_Sel_pData_merge.RData")
write.csv(Table_top2_Sel,file = "Table_top2_Sel.csv")

# For top2 (delta betas =>20%)
# Merge top_sva with EPIC_anno
de <- merge(top2, EPIC_anno, by="row.names",all.x= T)
write.table(de, file="top2_ComBat_Batch_annotated_050718.txt", row.names=T, sep="\t")



# Make heatmap for the selected sites
# read the heatmap making file
heatmap<-read.csv("Table_top2_Sel_For_Heatmap_Gene_CpG.csv", header=T)
heatmap2 <- heatmap[,-1]
rownames(heatmap2) <- heatmap[,1]
heatmap2<-as.matrix(heatmap2)
colnames(heatmap2)

tiff("Heatmap_ComBatcorrected.tif", res=300, compression = "lzw", height=30, width=30, units="in")
heatmap(heatmap2,cexRow =1,cexCol=3)
dev.off()

tiff("Heatmap_ComBatcorrected_2.tif", res=300, compression = "lzw", height=30, width=10, units="in")
heatmap(heatmap2, Rowv=NA ,cexRow =1,cexCol=3)
dev.off()

library(gplots)
#or try heatplot 
tiff("Heatmap_ComBatcorrected_3.tif", res=300, compression = "lzw", height=40, width=40, units="in")
heatmap.2(heatmap2, scale="row",dendrogram="col",cexRow =1,cexCol=3)
dev.off()

# Volcano plot
tiff("Volcano_Plot_3.tif", res=300, compression = "lzw", height=10, width=10, units="in")
volcanoplot(eb_sva,coef=1,highlight=row.names(top_sva))
dev.off()

############################################

# Pie charts to see Hyper and hypomethylated CpGs in top_sva (<=0.05)
# and top_sva (<=0.05, DeltaBetas_20)

############################################
# All significant CpGs
x <-  c(8379,3424)
labels <-  c("Hypomethylated","Hypermethylated")
piepercent<- round(100*x/sum(x), 1)
tiff("Hypo-Hyper_All_Pie Chart.tif", res=300, compression = "lzw", height=5, width=10, units="in")
pie(x, labels = piepercent, main = "Hypo-Hyper_All_pie chart",col = rainbow(length(x)))
legend("topright", c("Hypomethylated","Hypermethylated"), cex = 0.8,
       fill = rainbow(length(x)))
dev.off()

# All significant CpGs with Delta Betas => 20%
x <-  c(162,94)
labels <-  c("Hypomethylated","Hypermethylated")
piepercent<- round(100*x/sum(x), 1)
tiff("Hypo-Hyper_Delta_Beta20_Pie Chart.tif", res=300, compression = "lzw", height=5, width=10, units="in")
pie(x, labels = piepercent, main = "Hypo-Hyper_Delta_Beta20_Pie chart",col = rainbow(length(x)))
legend("topright", c("Hypomethylated","Hypermethylated"), cex = 0.8,
       fill = rainbow(length(x)))
dev.off()



############################################
# GO enrichment and pathway analysis using missMethyl
############################################
library(missMethyl)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
library(topGO)
load("mSetSqFlt_CR_280618.RData")
rownames(mSetSqFlt_CR)<-as.character(rownames(mSetSqFlt_CR))
# Load significant CpGs (11803 CpG sites)
top<-read.table("top_ComBat_Batch_050718.txt", header=T)
top<-as.character(rownames(top))

###################################
# GO analysis
##################################
gst <- gometh(sig.cpg=top, all.cpg=rownames(mSetSqFlt_CR),array.type= "EPIC", collection="GO")
newdata <- gst[order(gst$FDR),]

###################################
# KEGG pathway analysis
##################################
# Without probe bias
gst_KEGG <- gometh(sig.cpg=top, all.cpg=rownames(mSetSqFlt_CR),array.type= "EPIC", collection="KEGG")
newdata_KEGG_Allsignificant <- gst_KEGG[order(gst_KEGG$FDR),]
write.csv (newdata_KEGG_Allsignificant, "KEGG_pathways_Allsignificant_090718.csv")

############################################
# Bias Plot
############################################
pdf(file="Bias_Plot_All_significant.pdf")
gst <- gometh(sig.cpg=top, all.cpg=rownames(mSetSqFlt_CR),array.type= "EPIC",plot.bias=TRUE,prior.prob = TRUE)
dev.off()


#######################################
# For 256 CpGs with delta beta of=>20%
#######################################
top<-read.table("top2_ComBat_Batch_deltaBetas20_050718.txt", header=T)
top<-as.character(rownames(top))
gst <- gometh(sig.cpg=top, all.cpg=rownames(mSetSqFlt_CR),array.type= "EPIC", collection="GO")
newdata <- gst[order(gst$FDR),]

######################
# KEGG pothway analysis
######################
gst_KEGG <- gometh(sig.cpg=top, all.cpg=rownames(mSetSqFlt_CR),array.type= "EPIC", collection="KEGG")
newdata_KEGG <- gst_KEGG[order(gst_KEGG$FDR),]
write.csv (newdata_KEGG, "KEGG_pathways_DeltaBetas_090718.csv")

############################################
# Bias Plot
############################################
pdf(file="Bias_Plot_DeltaBetas_20.pdf")
gst <- gometh(sig.cpg=top, all.cpg=rownames(mSetSqFlt_CR),array.type= "EPIC",plot.bias=TRUE,prior.prob = TRUE)
dev.off()


####################################################################
# DNA methylation of Podocyte-Associated Molecules in response to IR
####################################################################

load("mSetSqFlt_CR_280618.RData")
meth <- getM(mSetSqFlt_CR)
betas<-ilogit2(meth)
rownames(betas)
slit_genes<-c("cg18374265","cg20231769","cg22339278", "cg24301250","cg01337940", "cg01719390", "cg13991728")
NPHS1_Meth<-betas[rownames(betas)%in% slit_genes,]
# Transpose and convert to data.frame
NPHS1_Meth<-as.data.frame(NPHS1_Meth)
t_NPHS1_Meth<-t(NPHS1_Meth)
t_NPHS1_Meth<-data.frame(t_NPHS1_Meth)
names(t_NPHS1_Meth)[1]<-paste("SYNPO2_Meth")
names(t_NPHS1_Meth)[2]<-paste("FAT1_Meth")
names(t_NPHS1_Meth)[3]<-paste("SYNPO_Meth")
names(t_NPHS1_Meth)[4]<-paste("KIRREL3_Meth")
names(t_NPHS1_Meth)[5]<-paste("NPHS1_Meth")
names(t_NPHS1_Meth)[6]<-paste("KIRREL2_A_Meth")
names(t_NPHS1_Meth)[7]<-paste("KIRREL2_B_Meth")

t_NPHS1_Meth$Barcode<-rownames(t_NPHS1_Meth)
# Load the phenotype data
target<-read.csv("Sample_sheet.csv", header=T)
# merge with phenotype Data
merge<-merge(t_NPHS1_Meth, target, by.x= "Barcode", by.y= "Basename")
# Boxplot 
library(reshape2)
# Make a dtaframe with only methylation values for NPHS1 and Group
test<-merge[,c("Group","SYNPO2_Meth","FAT1_Meth","SYNPO_Meth","KIRREL3_Meth","NPHS1_Meth","KIRREL2_A_Meth","KIRREL2_B_Meth")]
test.m <- melt(test)
library(ggplot2)
tiff("boxplot_SYNPO2_SYNPO_NPHS_KIRREL_FAT1_Methylation_IR_IS_Podo.tif", res=300, compression = "lzw", height=5, width=15, units="in")
ggplot(test.m, aes(x=variable, y=value, fill =  Group)) +
  geom_boxplot(outlier.colour = "black", outlier.size = NA) +
  scale_y_continuous(limits = quantile(test.m$value, c(0.01, 1.00))) +
  scale_fill_manual(values = c("coral2", "cornflowerblue"))+
  xlab("Insulin Resistance Response") +
  ylab("DNA methylation (Beta values)")+
  ggtitle("")
dev.off()  

####################################################################
# DNA methylation of WT1 and Beta catenin in response to IR
####################################################################

load("mSetSqFlt_CR_280618.RData")
meth <- getM(mSetSqFlt_CR)
betas<-ilogit2(meth)
rownames(betas)
slit_genes<-c("cg00142072","cg25563456","cg20914573")
WT1_Meth<-betas[rownames(betas)%in% slit_genes,]
# Transpose and convert to data.frame
WT1_Meth<-as.data.frame(WT1_Meth)
t_WT1_Meth<-t(WT1_Meth)
t_WT1_Meth<-data.frame(t_WT1_Meth)
names(t_WT1_Meth)[1]<-paste("CTNNB1_Meth_Body")
names(t_WT1_Meth)[2]<-paste("WT1_Meth_TSS200")
names(t_WT1_Meth)[3]<-paste("WT1_Meth_Body")


t_WT1_Meth$Barcode<-rownames(t_WT1_Meth)
# Load the phenotype data
target<-read.csv("Sample_sheet.csv", header=T)
# merge with phenotype Data
merge<-merge(t_WT1_Meth, target, by.x= "Barcode", by.y= "Basename")
# Boxplot 
library(reshape2)
# Make a dataframe with only methylation values for NPHS1 and Group
test<-merge[,c("Group","CTNNB1_Meth_Body","WT1_Meth_TSS200","WT1_Meth_Body")]
test.m <- melt(test)
library(ggplot2)
tiff("boxplot_WT1_CTNNB1_Methylation_IR_IS_Podo.tif", res=300, compression = "lzw", height=5, width=15, units="in")
ggplot(test.m, aes(x=variable, y=value, fill =  Group)) +
  geom_boxplot(outlier.colour = "black", outlier.size = NA) +
  scale_y_continuous(limits = quantile(test.m$value, c(0.01, 1.00))) +
  scale_fill_manual(values = c("coral2", "cornflowerblue"))+
  xlab("Insulin Resistance Response") +
  ylab("DNA methylation (Beta values)")+
  ggtitle("")
dev.off()  


####################################################################
# DNA methylation of ABCC4, MET and SLC2A5 in response to IR (From paper by Takeshi Marumo et al. JASN 2015)
# They used mouse model
####################################################################
load("mSetSqFlt_CR_280618.RData")
meth <- getM(mSetSqFlt_CR)
betas<-ilogit2(meth)
rownames(betas)
slit_genes<-c("cg24007802","cg08964563","cg17341110")
WT1_Meth<-betas[rownames(betas)%in% slit_genes,]
# Transpose and convert to data.frame
WT1_Meth<-as.data.frame(WT1_Meth)
t_WT1_Meth<-t(WT1_Meth)
t_WT1_Meth<-data.frame(t_WT1_Meth)
names(t_WT1_Meth)[1]<-paste("SLC2A5_Meth_Body")
names(t_WT1_Meth)[2]<-paste("MET_Meth_5UTR")
names(t_WT1_Meth)[3]<-paste("ABCC4_Meth_Body")


t_WT1_Meth$Barcode<-rownames(t_WT1_Meth)
# Load the phenotype data
target<-read.csv("Sample_sheet.csv", header=T)
# merge with phenotype Data
merge<-merge(t_WT1_Meth, target, by.x= "Barcode", by.y= "Basename")
# Boxplot 
library(reshape2)
# Make a dataframe with only methylation values for NPHS1 and Group
test<-merge[,c("Group","SLC2A5_Meth_Body","MET_Meth_5UTR","ABCC4_Meth_Body")]
test.m <- melt(test)
library(ggplot2)
tiff("boxplot_SLC2A5_MET_ABCC4_Methylation_IR_IS_Podo.tif", res=300, compression = "lzw", height=5, width=10, units="in")
ggplot(test.m, aes(x=variable, y=value, fill =  Group)) +
  geom_boxplot(outlier.colour = "black", outlier.size = NA) +
  scale_y_continuous(limits = quantile(test.m$value, c(0.01, 1.00))) +
  scale_fill_manual(values = c("coral2", "cornflowerblue"))+
  xlab("Insulin Resistance Response") +
  ylab("DNA methylation (Beta values)")+
  ggtitle("")
dev.off()  



####################################################################
# DNA methylation of Top 10 CpG in response to IR (Delta Betas =>20%, FDR <=0.05)
#
####################################################################
load("mSetSqFlt_CR_280618.RData")
meth <- getM(mSetSqFlt_CR)
betas<-ilogit2(meth)
rownames(betas)
top_CpG<-c("cg00945023","cg01427750","cg02283643","cg02366189","cg02403137",
           "cg02438545","cg03025760","cg03177023","cg03254818","cg03290514")
top_Meth<-betas[rownames(betas)%in% top_CpG,]
# Transpose and convert to data.frame
top_Meth<-as.data.frame(top_Meth)
t_top_Meth<-t(top_Meth)
t_top_Meth<-data.frame(t_top_Meth)
names(t_top_Meth)[1]<-paste("C2orf80_Meth_TSS200")
names(t_top_Meth)[2]<-paste("SRD5A3-AS1_Meth_TSS200")
names(t_top_Meth)[3]<-paste("cg03177023_Meth_Intergenic")
names(t_top_Meth)[4]<-paste("cg03254818_Meth_Intergenic")
names(t_top_Meth)[5]<-paste("cg02403137_Meth_Intergenic")
names(t_top_Meth)[6]<-paste("SULF1_Meth_TSS200")
names(t_top_Meth)[7]<-paste("CDC14B_Meth_Body")
names(t_top_Meth)[8]<-paste("GPR133_Meth_Body")
names(t_top_Meth)[9]<-paste("cg02366189_Meth_Intergenic")
names(t_top_Meth)[10]<-paste("SLITRK2_Meth_3'UTR")


t_top_Meth$Barcode<-rownames(t_top_Meth)
# Load the phenotype data
target<-read.csv("Sample_sheet.csv", header=T)
# merge with phenotype Data
merge<-merge(t_top_Meth, target, by.x= "Barcode", by.y= "Basename")
# Boxplot 
library(reshape2)
# Make a dataframe with only methylation values for NPHS1 and Group
test<-merge[,c("Group","C2orf80_Meth_TSS200","SRD5A3-AS1_Meth_TSS200","cg03177023_Meth_Intergenic",
               "cg03254818_Meth_Intergenic", "cg02403137_Meth_Intergenic", "SULF1_Meth_TSS200",
               "CDC14B_Meth_Body","GPR133_Meth_Body","cg02366189_Meth_Intergenic", "SLITRK2_Meth_3'UTR" )]
test.m <- melt(test)
library(ggplot2)
tiff("boxplot_Top_10_Methylation_IR_IS_Podo.tif", res=300, compression = "lzw", height=10, width=20, units="in")
ggplot(test.m, aes(x=variable, y=value, fill =  Group)) +
  geom_boxplot(outlier.colour = "black", outlier.size = NA) +
  scale_y_continuous(limits = quantile(test.m$value, c(0.01, 1.00))) +
  scale_fill_manual(values = c("coral2", "cornflowerblue"))+
  xlab("Insulin Resistance Response") +
  ylab("DNA methylation (Beta values)")+
  ggtitle("")
dev.off() 





# Load the EPIC_Annotation for EPIC
rownames(EPIC_anno)<-EPIC_anno$name
Beta_Anno<-EPIC_anno[rownames(EPIC_anno)%in% rownames(meth),]
colnames(Beta_Anno)
head(Beta_Anno)
Beta_Anno_Reg<-Beta_Anno[,c("gene.region","relation.to.island")]
write.csv(Beta_Anno_Reg, file="Beta_Anno_Reg.csv")


# Making Pie-chart for Significant CpGs, Delta betas (>20%) and Overall- For Gene centric probes

# For All significant CpGs
x <-  c(303,329,952,4619,1833,3767)
labels <-  c("1stExon","3'UTR","5'UTR","Body","TSS","Intergenic")
piepercent<- round(100*x/sum(x), 1)
tiff("Gene Centric Pie Chart.tif", res=300, compression = "lzw", height=5, width=10, units="in")
pie(x, labels = piepercent, main = "Gene Centric pie chart",col = rainbow(length(x)))
legend("topright", c("1stExon","3'UTR","5'UTR","Body","TSS","Intergenic"), cex = 0.8,
       fill = rainbow(length(x)))
dev.off()

# For Delta betas more than 20% CpGs
x <-  c(7,10,17,98,31,93)
labels <-  c("1stExon","3'UTR","5'UTR","Body","TSS","Intergenic")
piepercent<- round(100*x/sum(x), 1)
tiff("Gene Centric Pie Chart Delta Betas 20.tif", res=300, compression = "lzw", height=5, width=10, units="in")
pie(x, labels = piepercent, main = "Gene Centric pie chart",col = rainbow(length(x)))
legend("topright", c("1stExon","3'UTR","5'UTR","Body","TSS","Intergenic"), cex = 0.8,
       fill = rainbow(length(x)))
dev.off()

# For All CpGs
x <-  c(29694,19232,66916,283931,158983,213734)
labels <-  c("1stExon","3'UTR","5'UTR","Body","TSS","Intergenic")
piepercent<- round(100*x/sum(x), 1)
tiff("Gene Centric Pie Chart_850k.tif", res=300, compression = "lzw", height=5, width=10, units="in")
pie(x, labels = piepercent, main = "Gene Centric pie chart",col = rainbow(length(x)))
legend("topright", c("1stExon","3'UTR","5'UTR","Body","TSS","Intergenic"), cex = 0.8,
       fill = rainbow(length(x)))
dev.off()


# Making Pie-chart for Significant CpGs, Delta betas (>20%) and Overall- For CpG Island probes

# For All significant CpGs
x <-  c(8107,750,1728,1218)
labels <-  c("Open Sea","Shelf","Shore","Island")
piepercent<- round(100*x/sum(x), 1)
tiff("CpG Island Pie Chart.tif", res=300, compression = "lzw", height=5, width=10, units="in")
pie(x, labels = piepercent, main = "Gene Centric pie chart",col = rainbow(length(x)))
legend("topright", c("Open Sea","Shelf","Shore","Island"), cex = 0.8,
       fill = rainbow(length(x)))
dev.off()

# For Delta betas more than 20% CpGs
x <-  c(146,21,47,42)
labels <-  c("Open Sea","Shelf","Shore","Island")
piepercent<- round(100*x/sum(x), 1)
tiff("CpG Island Pie Chart Delta Betas 20.tif", res=300, compression = "lzw", height=5, width=10, units="in")
pie(x, labels = piepercent, main = "Gene Centric pie chart",col = rainbow(length(x)))
legend("topright", c("Open Sea","Shelf","Shore","Island"), cex = 0.8,
       fill = rainbow(length(x)))
dev.off()

# For All CpGs
x <-  c(425721,53481,141659,151629)
labels <-  c("Open Sea","Shelf","Shore","Island")
piepercent<- round(100*x/sum(x), 1)
tiff("CpG Island Pie Chart_850k.tif", res=300, compression = "lzw", height=5, width=10, units="in")
pie(x, labels = piepercent, main = "Gene Centric pie chart",col = rainbow(length(x)))
legend("topright", c("Open Sea","Shelf","Shore","Island"), cex = 0.8,
       fill = rainbow(length(x)))
dev.off()



################################################
# DMR analysis using Ewaff
# Done on Bluecrystal
################################################
install.packages("/panfs/panasas01/sscm/sa16666/450kanalysis/Podocyte_Methylation/ewaff.tar.gz", repos = NULL, type="source")
head(eb_sva$coefficients)
estimate<-eb_sva$coefficients
estimate<-as.data.frame(estimate)
#estimate<-estimate[,1]
se<-eb_sva$stdev.unscaled * eb_sva$sigma
head(se)
se<-as.data.frame(se)
#se<-se[,1]
p.value<-eb_sva$p.value
head(p.value)
#p.value<-p.value[,1]
p.value<-as.data.frame(p.value)

#merge three data frames
df <- cbind.data.frame(estimate,se,p.value)

names(df)[1]<-paste("estimate")
names(df)[2]<-paste("se")
names(df)[3]<-paste("p.value")

load("EPIC_anno_meffil.RData")
rownames(EPIC_anno)<-EPIC_anno$name
Beta_Anno<-EPIC_anno[rownames(EPIC_anno)%in% rownames(meth),]

Beta_Anno<-Beta_Anno[,c("name","chromosome","position")]
Beta_Anno_2<-Beta_Anno[order(Beta_Anno$chromosome,Beta_Anno$position),]

df2<-df[rownames(df) %in% rownames(Beta_Anno),]

combat_meth_2<-combat_meth[rownames(combat_meth) %in% rownames(Beta_Anno),]


bumps<-ewaff.bumps(estimate=df2$estimate, se=df2$se, p.value=df2$p.value, methylation=combat_meth_2, chr=Beta_Anno_2$chromosome , pos=Beta_Anno_2$position, maxgap = 500,
            p.cutoff = 0.05, verbose = T)

bumps_FDR<-bumps[which(bumps$p.adjust < 0.05),
      c("chr","start","end","n", "p.value","p.adjust")]


write.csv(bumps_FDR, file="bumps_FDR_Podocytes.csv")


#######################################################
# DMR analysis using ChAMP
######################################################
source("https://bioconductor.org/biocLite.R")
biocLite("ChAMP")
library(ChAMP)
library(sva)

load("mSetSqFlt_CR_280618.RData")
meth <- getM(mSetSqFlt_CR)
head(meth[,1:5])
meth[!is.finite(meth)] <- 0
# Load phenotype data
target<-read.csv("Sample_sheet.csv", header=T)


# ComBat correct methylation data
batch = target$Batch
combat_meth = ComBat(dat=meth, batch=batch, mod=NULL, par.prior=TRUE, prior.plots=FALSE)
#combat_meth<-as.matrix(combat_meth)
#mod<-model.matrix(~as.factor(Group), data=target)

load("EPIC_anno_meffil.RData")
head(EPIC_anno)
Beta_Anno<-EPIC_anno[,c(3,4,5)]
rownames(Beta_Anno)<-Beta_Anno$name

rownames(combat_meth)
rownames(Beta_Anno)

# Merge Beta Anno and Meth
combat_meth<-as.data.frame(combat_meth)
Beta_Anno<-as.data.frame(Beta_Anno)

Merge<-merge(combat_meth,Beta_Anno,by="row.names",all.x=TRUE)

rownames(EPIC_anno)<-EPIC_anno$name
Beta_Anno<-EPIC_anno[rownames(EPIC_anno)%in% rownames(meth),]
colnames(Beta_Anno)
rownames(Beta_Anno)
# convert combat_meth into a matrix
combat_meth<-as.matrix(combat_meth)

myDMR <- champ.DMR(beta=combat_meth,pheno=target$Group,method="Bumphunter",arraytype= "EPIC",nullMethod="bootstrap",B=250)

head(myDMR$BumphunterDMR)

write.csv(myDMR, file="Bumphunter_Podocytes.csv")



###########################################
# Model wtihout ComBat correction (None of the CpG sites reached statistical significance)
###########################################
# See possible effects of technical biological factors before running the analysis:
mod<-model.matrix(~as.factor(Group), data=target)
fit_sva= lmFit(meth,mod)
contrast.matrix<-c(0,1)
fitContrasts_sva = contrasts.fit(fit_sva,contrast.matrix)
eb_sva = eBayes(fitContrasts_sva)
top_sva<-topTable(eb_sva, adjust="BH",number=1000000,p.value=0.3,genelist=eb_sva$genes)
# Write table with FDR <=0.05



###################################################

# MR analysis of our findings

####################################################
# Load the mQTLs from ARIES database for 11803 CpGs
# Done manually
# Load the data for mQTL

setwd("C:/MY documents/Podocyte_methylation")
newdata<-read.csv("SNP_CpG_11803_Podocytes.csv", header=T)

# Calculate SE from p value and beta
newdata$se=abs(newdata$beta/qnorm(newdata$p.value/2))

# To confirm if we have derived the right SE from the dataset, we rederive the P values from SE as follows:
newdata$pvalue_Re=2*pnorm(abs(newdata$beta/newdata$se), lower.tail=F)
# It works!

# Remove Trans SNPs

newdata<-newdata[!(newdata$Trans=="Y"),]


# Remove columns not needed for further analysis and save as CSV file 
newdata_1 <- newdata[ -c(1,7,8,9,10,12,13,15,17) ]

names(newdata_1)[2]<-"CHR"
names(newdata_1)[3]<-"SNP_Pos"
names(newdata_1)[4]<-"effect_allele"
names(newdata_1)[5]<-"other_allele"
names(newdata_1)[7]<-"pval"



levels(newdata_1$effect_allele)

# Remove the rows if effect allele is "D", "I","R"
# 
newdata_2<-newdata_1[!(newdata_1$effect_allele=="D"),]

newdata_3<-newdata_2[!(newdata_2$effect_allele=="I"),]

newdata_4<-newdata_3[!(newdata_3$effect_allele=="R"),]


write.table(newdata_4, "ARIES_Childhood_11803_CpGs_Cis_for_MR.txt", sep="\t")


# Perfom MR analysis
# Carry out clumping
library(TwoSampleMR)

file <- read.table("For_Clumping_ARIES_CHR_1_CpGs_Cis_for_MR.txt", header=T)
head(file)



file_cl_1 <- clump_data(file,clump_r2 = 0.001)


file <- read.table("For_Clumping_ARIES_CHR_2_CpGs_Cis_for_MR.txt", header=T)
head(file)
file_cl_2 <- clump_data(file,clump_r2 = 0.001)


file <- read.table("For_Clumping_ARIES_CHR_3_CpGs_Cis_for_MR.txt", header=T)
head(file)
file_cl_3 <- clump_data(file,clump_r2 = 0.001)


file <- read.table("For_Clumping_ARIES_CHR_4_CpGs_Cis_for_MR.txt", header=T)
head(file)
file_cl_4 <- clump_data(file,clump_r2 = 0.001)


file <- read.table("For_Clumping_ARIES_CHR_5_CpGs_Cis_for_MR.txt", header=T)
head(file)
file_cl_5 <- clump_data(file,clump_r2 = 0.001)

file <- read.table("For_Clumping_ARIES_CHR_6_CpGs_Cis_for_MR.txt", header=T)
head(file)
file_cl_6 <- clump_data(file,clump_r2 = 0.001)


file <- read.table("For_Clumping_ARIES_CHR_7_CpGs_Cis_for_MR.txt", header=T)
head(file)
file_cl_7 <- clump_data(file,clump_r2 = 0.001)


file <- read.table("For_Clumping_ARIES_CHR_8_CpGs_Cis_for_MR.txt", header=T)
head(file)
file_cl_8 <- clump_data(file,clump_r2 = 0.001)

file <- read.table("For_Clumping_ARIES_CHR_9_CpGs_Cis_for_MR.txt", header=T)
head(file)
file_cl_9 <- clump_data(file,clump_r2 = 0.001)

file <- read.table("For_Clumping_ARIES_CHR_10_CpGs_Cis_for_MR.txt", header=T)
head(file)
file_cl_10 <- clump_data(file,clump_r2 = 0.001)


file <- read.table("For_Clumping_ARIES_CHR_11_CpGs_Cis_for_MR.txt", header=T)
head(file)
file_cl_11 <- clump_data(file,clump_r2 = 0.001)


file <- read.table("For_Clumping_ARIES_CHR_12_CpGs_Cis_for_MR.txt", header=T)
head(file)
file_cl_12 <- clump_data(file,clump_r2 = 0.001)

file <- read.table("For_Clumping_ARIES_CHR_13_CpGs_Cis_for_MR.txt", header=T)
head(file)
file_cl_13 <- clump_data(file,clump_r2 = 0.001)


file <- read.table("For_Clumping_ARIES_CHR_14_CpGs_Cis_for_MR.txt", header=T)
head(file)
file_cl_14 <- clump_data(file,clump_r2 = 0.001)


file <- read.table("For_Clumping_ARIES_CHR_15_CpGs_Cis_for_MR.txt", header=T)
head(file)
file_cl_15 <- clump_data(file,clump_r2 = 0.001)


file <- read.table("For_Clumping_ARIES_CHR_16_CpGs_Cis_for_MR.txt", header=T)
head(file)
file_cl_16 <- clump_data(file,clump_r2 = 0.001)


file <- read.table("For_Clumping_ARIES_CHR_17_CpGs_Cis_for_MR.txt", header=T)
head(file)
file_cl_17 <- clump_data(file,clump_r2 = 0.001)

file <- read.table("For_Clumping_ARIES_CHR_18_CpGs_Cis_for_MR.txt", header=T)
head(file)
file_cl_18 <- clump_data(file,clump_r2 = 0.001)

file <- read.table("For_Clumping_ARIES_CHR_19_CpGs_Cis_for_MR.txt", header=T)
head(file)
file_cl_19 <- clump_data(file,clump_r2 = 0.001)

file <- read.table("For_Clumping_ARIES_CHR_20_CpGs_Cis_for_MR.txt", header=T)
head(file)
file_cl_20 <- clump_data(file,clump_r2 = 0.001)


file <- read.table("For_Clumping_ARIES_CHR_21_CpGs_Cis_for_MR.txt", header=T)
head(file)
file_cl_21 <- clump_data(file,clump_r2 = 0.001)


file <- read.table("For_Clumping_ARIES_CHR_22_CpGs_Cis_for_MR.txt", header=T)
head(file)
file_cl_22 <- clump_data(file,clump_r2 = 0.001)



## Row bind all the dataframes

library(gtools)
MR_SNP<-smartbind(file_cl_22,file_cl_21,file_cl_20,file_cl_19,file_cl_18,file_cl_17,file_cl_16,file_cl_15,file_cl_14,file_cl_13,file_cl_12,file_cl_11,file_cl_10,file_cl_9,file_cl_8,file_cl_7,file_cl_6,file_cl_5,file_cl_4,file_cl_3,file_cl_2,file_cl_1)

write.csv(MR_SNP, file="MR_SNP_clumped_Podocytes.csv")

##########################
# Carry out the 2S MR
##########################

MR_SNP<-read.csv("MR_SNP_clumped_Podocytes.csv", header=T)


exposure_dat<-read_exposure_data(filename="MR_SNP_clumped_Podocytes.csv", clump = FALSE, sep = ",",
                   phenotype_col = "Phenotype", snp_col = "SNP", beta_col = "beta",
                   se_col = "se", effect_allele_col = "effect_allele",
                   other_allele_col = "other_allele", pval_col = "pval", min_pval = 1e-500)



ao <- available_outcomes()

subset_ao<-subset(ao, select=c(trait, id))

# Get effects of instruments on outcome
# For Chronic Kidney Disease
outcome_dat <- extract_outcome_data(snps=exposure_dat$SNP, outcomes=1102)
# Harmonise the exposure and outcome data
dat <- harmonise_data(exposure_dat, outcome_dat)

# Perform MR
res <- mr(dat)

res

write.csv(res, file="Results_MR_CKD_mQTL_Podocytes.csv")


# For Urinary Albumin to Creatinine Ratio
outcome_dat <- extract_outcome_data(snps=exposure_dat$SNP, outcomes=1100)

# Harmonise the exposure and outcome data
dat <- harmonise_data(exposure_dat, outcome_dat)

# Perform MR
res <- mr(dat)

res

write.csv(res, file="Results_MR_Alb_Crea_Ratio_mQTL_Podocytes.csv")

# For Urinary Albumin to Creatinine Ratio
outcome_dat <- extract_outcome_data(snps=exposure_dat$SNP, outcomes=1101)

# Harmonise the exposure and outcome data
dat <- harmonise_data(exposure_dat, outcome_dat)

# Perform MR
res <- mr(dat)

res

write.csv(res, file="Results_II_MR_Alb_Crea_Ratio_mQTL_Podocytes.csv")


# For urinary Albumin to creatinine ratio (Not working !!)
outcome_dat <- extract_outcome_data(snps=exposure_dat$SNP, outcomes=1107)

# Harmonise the exposure and outcome data
dat <- harmonise_data(exposure_dat, outcome_dat)

# Perform MR
res <- mr(dat)

res

write.csv(res, file="Results_III_MR_Alb_Crea_Ratio_mQTL_Podocytes.csv")


# For Microalbuminuria
outcome_dat <- extract_outcome_data(snps=exposure_dat$SNP, outcomes=1097)

# Harmonise the exposure and outcome data
dat <- harmonise_data(exposure_dat, outcome_dat)

# Perform MR
res <- mr(dat)

res

write.csv(res, file="Results_MR_Microalbuminuria_mQTL_Podocytes.csv")



# For serum creatinine eGFR
outcome_dat <- extract_outcome_data(snps=exposure_dat$SNP, outcomes=1103)

# Harmonise the exposure and outcome data
dat <- harmonise_data(exposure_dat, outcome_dat)

# Perform MR
res <- mr(dat)

res

write.csv(res, file="Results_I_MR_Serum creatinine_eGFR_mQTL_Podocytes.csv")


outcome_dat <- extract_outcome_data(snps=exposure_dat$SNP, outcomes=1104)

# Harmonise the exposure and outcome data
dat <- harmonise_data(exposure_dat, outcome_dat)

# Perform MR
res <- mr(dat)

res

write.csv(res, file="Results_II_MR_Serum creatinine_eGFR_mQTL_Podocytes.csv")

outcome_dat <- extract_outcome_data(snps=exposure_dat$SNP, outcomes=1105)

# Harmonise the exposure and outcome data
dat <- harmonise_data(exposure_dat, outcome_dat)

# Perform MR
res <- mr(dat)

res

write.csv(res, file="Results_III_MR_Serum creatinine_eGFR_mQTL_Podocytes.csv")


# For T2D
outcome_dat <- extract_outcome_data(snps=exposure_dat$SNP, outcomes=1090)

# Harmonise the exposure and outcome data
dat <- harmonise_data(exposure_dat, outcome_dat)

# Perform MR
res <- mr(dat)

res

write.csv(res, file="Results_MR_T2D_mQTL_Podocytes.csv")

# For T2D
outcome_dat <- extract_outcome_data(snps=exposure_dat$SNP, outcomes=23)

# Harmonise the exposure and outcome data
dat <- harmonise_data(exposure_dat, outcome_dat)

# Perform MR
res <- mr(dat)

res

write.csv(res, file="Results_II_MR_T2D_mQTL_Podocytes.csv")



# Fasting insulin
outcome_dat <- extract_outcome_data(snps=exposure_dat$SNP, outcomes=775)

# Harmonise the exposure and outcome data
dat <- harmonise_data(exposure_dat, outcome_dat)

# Perform MR
res <- mr(dat)

res

write.csv(res, file="Results_Fasting_Insulin_mQTL_Podocytes.csv")


# For HOMA-IR
outcome_dat <- extract_outcome_data(snps=exposure_dat$SNP, outcomes=771)

# Harmonise the exposure and outcome data
dat <- harmonise_data(exposure_dat, outcome_dat)

# Perform MR
res <- mr(dat)

res

write.csv(res, file="Results_HOMA-IR_mQTL_Podocytes.csv")



# For Coronary Heart Disease (Not working !!)
outcome_dat <- extract_outcome_data(snps=exposure_dat$SNP, outcomes=9)
# Harmonise the exposure and outcome data
dat <- harmonise_data(exposure_dat, outcome_dat)

# Perform MR
res <- mr(dat)

res

write.csv(res, file="Results_CHD_mQTL_Podocytes.csv")

