### Load Libraries
##################
library(chipseq)
library(GenomicRanges)
library(htSeqTools)
library(rtracklayer)
library(spp)
library(limma)
library(ChIPpeakAnno)
##################
### Get data
##################
ChIPDirectory <- "/Data_For_ChIP_Practical/"
BamFiles <- dir(file.path(getwd(),ChIPDirectory),pattern="*.bam$",full.name=T)

##Make the names a little tidier
names(BamFiles) <- gsub("Chr1|\\.bam","",basename(BamFiles))
##Put into a BamFileList for use later on.
BamsForAnalysis <- BamFileList(BamFiles)
##################
### Access Header
##################
print(path(BamsForAnalysis))
samHeader <- scanBamHeader(path(BamsForAnalysis["TF_1"]))
str(samHeader,max.level=1)
print(samHeader[[1]]$targets)
print(samHeader[[1]]$text)
# Sort Order
print(samHeader[[1]]$text["@HD"])
# Program used in alignment
print(samHeader[[1]]$text["@PG"])
# Species
print(samHeader[[1]]$text["@CO"])

AllSeqnames <- samHeader[[1]]$targets
LongestChr <- AllSeqnames[names(AllSeqnames) == "chr1"]
LCName <- names(LongestChr)
LCLength <- unique(LongestChr)

?scanBamFlag

which <- GRanges(seqnames=LCName,ranges=IRanges(1,LCLength))
param <- ScanBamParam(which=which)
param
##################
### Read in the Data
##################
BamBnd <-  readGappedAlignments(path(BamsForAnalysis["TF_1"]),param=param)
GrangesAlign <- granges(BamBnd)
seqlevels(GrangesAlign) <- "chr1"
ReadLength <- round(median(width(GrangesAlign[1:1000])))

GrangesAlign[1:5]

##################
### Make Coverage plots
##################

AlignReadEnds_Pos <- GrangesAlign[strand(GrangesAlign) == "+"]
AlignReadEnds_Neg <- GrangesAlign[strand(GrangesAlign) == "-"]

All_Coverage <- coverage(GrangesAlign)
Pos_Coverage <- coverage(AlignReadEnds_Pos)
Neg_Coverage <- coverage(AlignReadEnds_Neg)

print(All_Coverage)

Range <- c(211646604:211649812)

Graph_RegionOfInterest_Pos <- Pos_Coverage$chr1[211646604:211649812]
Graph_RegionOfInterest_Neg <- Neg_Coverage$chr1[211646604:211649812]
smoothedNegCoverage <- caTools::runmean(as.vector(Graph_RegionOfInterest_Neg),20)
smoothedPosCoverage <- caTools::runmean(as.vector(Graph_RegionOfInterest_Pos),20)

png("NoReadExtension.png")
plot(smoothedNegCoverage,col="red",type="l",ylab="Smoothed Read Depth 20bp Window",xlab="bp from range start")
lines(smoothedPosCoverage,col="green",type="l")
dev.off()
which.max(smoothedNegCoverage) - which.max(smoothedPosCoverage)

##################
### Calculate Cross-Coverage and remake plots
##################


fraglenCov <- median(estimate.mean.fraglen(c(AlignReadEnds_Pos,AlignReadEnds_Neg),method="coverage"))
ExtendedReads <- resize(GrangesAlign,fraglenCov,fix="start")
ExtendedReads_Pos <- ExtendedReads[strand(ExtendedReads) == "+"]
ExtendedReads_Neg <- ExtendedReads[strand(ExtendedReads) == "-"]


Graph_RegionOfInterest_Pos <- coverage(ExtendedReads_Pos)$chr1[Range]
Graph_RegionOfInterest_Neg <- coverage(ExtendedReads_Neg)$chr1[Range]
smoothedNegCoverage <- caTools::runmean(as.vector(Graph_RegionOfInterest_Neg),20)
smoothedPosCoverage <- caTools::runmean(as.vector(Graph_RegionOfInterest_Pos),20)

png("Cross_Coverage_ReadExtension.png")
plot(smoothedNegCoverage,col="red",type="l",ylab="Smoothed Read Depth 20bp Window",xlab="bp from range start")
lines(smoothedPosCoverage,col="green",type="l")
dev.off()
which.max(smoothedNegCoverage) - which.max(smoothedPosCoverage)

#################
### Cross Correlation and NSC/RSC
#################


library(spp)
readLength <- ReadLength

chip.data <-  read.bam.tags(path(BamsForAnalysis)[2])
binding.characteristics <- get.binding.characteristics(chip.data,srange=c(0,400),bin=2,accept.all.tags=T,remove.tag.anomalies = F)
crosscorr <- binding.characteristics


print(crosscorr$cross.correlation)



cc <- crosscorr$cross.correlation
MinY <- which.min(crosscorr$cross.correlation$y)
crosscorr$min.cc <- crosscorr$cross.correlation[MinY, ]
cat("Minimum cross-correlation value", crosscorr$min.cc$y,"\n",file=stdout())
cat("Minimum cross-correlation shift", crosscorr$min.cc$x,"\n",file=stdout())
sbw <- 2*floor(ceiling(5/2) / 2) + 1
cc$y <- runmean(cc$y,sbw,alg="fast")

bw <- ceiling(2/2)
peakidx <- (diff(cc$y,bw)>=0)
peakidx <- diff(peakidx,bw)
peakidx <- which(peakidx==-1) + bw

ShiftRangeExcludeMin <- 10
ShiftRangeExcludeMax <- readLength+10
peakidx <- peakidx[(cc$x[peakidx] < ShiftRangeExcludeMin) | (cc$x[peakidx] > ShiftRangeExcludeMax)]
cc <- cc[peakidx,]


maxpeakidx <- which.max(cc$y)
maxpeakshift <- cc$x[maxpeakidx]
maxpeakval <- cc$y[maxpeakidx]
peakidx <-which((cc$y >= 0.9*maxpeakval) & (cc$x >= maxpeakshift))
cc <- cc[peakidx,]

sortidx <- order(cc$y,decreasing=TRUE)
sortidx <- sortidx[c(1:min(3,length(sortidx)))]
cc.peak <- cc[sortidx,]

cat("Peak strand shift",paste(cc.peak$x,collapse=","),"\n",file=stdout())

crosscorr$peak$x <- cc.peak$x[1]
crosscorr$peak$y <- cc.peak$y[1]

ph.peakidx <- which( ( crosscorr$cross.correlation$x >= ( readLength - round(2*10) ) ) &
                     ( crosscorr$cross.correlation$x <= ( readLength + round(1.5*10) ) ) )
ph.peakidx <- ph.peakidx[ which.max(crosscorr$cross.correlation$y[ph.peakidx]) ]
crosscorr$phantom.cc <- crosscorr$cross.correlation[ph.peakidx,]
cat("Phantom peak location",crosscorr$phantom.cc$x,"\n",file=stdout())
cat("Phantom peak Correlation",crosscorr$phantom.cc$y,"\n",file=stdout())
crosscorr$phantom.coeff <- crosscorr$peak$y / crosscorr$phantom.cc$y
crosscorr$phantom.coeff <- crosscorr$peak$y / crosscorr$min.cc$y
cat("Normalized cross-correlation coefficient (NCCC)",crosscorr$phantom.coeff,"\n",file=stdout())
crosscorr$rel.phantom.coeff <- (crosscorr$peak$y - crosscorr$min.cc$y) / (crosscorr$phantom.cc$y - crosscorr$min.cc$y)
cat("Relative Cross correlation Coefficient (RCCC)",crosscorr$rel.phantom.coeff,"\n",file=stdout())
crosscorr$phantom.quality.tag <- NA

###################
## Now redo plot
###################

fraglenCor <- crosscorr$peak$x

ExtendedReads <- resize(GrangesAlign,fraglenCor,fix="start")

ExtendedReads_Pos <- ExtendedReads[strand(ExtendedReads) == "+"]
##Note the strand dependent start!!!
ExtendedReads_Neg <- ExtendedReads[strand(ExtendedReads) == "-"]


Graph_RegionOfInterest_Pos <- coverage(ExtendedReads_Pos)$chr1[Range]
Graph_RegionOfInterest_Neg <- coverage(ExtendedReads_Neg)$chr1[Range]
smoothedNegCoverage <- caTools::runmean(as.vector(Graph_RegionOfInterest_Neg),20)
smoothedPosCoverage <- caTools::runmean(as.vector(Graph_RegionOfInterest_Pos),20)

which.max(smoothedNegCoverage) - which.max(smoothedPosCoverage)
####################
## SSD/Gini and filtering Duplicates
####################

ssdOfSample <- ssdCoverage(GrangesAlign)
giniOfSample <- giniCoverage(GrangesAlign)

NumberOfDups <- tabDuplReads(RangedData(GrangesAlign))
CutOffs <- fdrEnrichedCounts(NumberOfDups,use=1:10,components=0,mc.cores=1)
DupFilt <- filterDuplReads(RangedData(GrangesAlign))
NumberOfDupsAfterEnrichmentTest <- tabDuplReads(DupFilt)
###############################################################

#####################
## Reading in peaks
#####################

TestMacsFile <- read.delim(file.path(getwd(),"/Data_For_ChIP_Practical/TF_1_peaks.bed"),sep="\t",header=F)
print(TestMacsFile[1:3,])

TestBed <- GRanges(seqnames=as.vector(TestMacsFile[,1]),IRanges(start=as.numeric(as.vector(TestMacsFile[,2])),end=as.numeric(as.vector(TestMacsFile[,3]))),strand=rep("*",nrow(TestMacsFile)))
elementMetadata(TestBed) <- TestMacsFile[,-c(1:3)]
colnames(elementMetadata(TestBed)) <- c("Peak_ID","Score")

TestBed[1:10]

length(TestBed)
width(TestBed)


##################
## Some useful  functions
##################
Bed2GRanges <- function(BedFileName){
    BedFile <- read.delim(BedFileName,sep="\t",header=F)
      StartPos <- 2
      EndPos <- 3
      ChrPos <- 1
      TempRanges_Bed <- GRanges(seqnames=as.vector(BedFile[,ChrPos]),IRanges(start=as.numeric(as.vector(BedFile[,StartPos])),end=as.numeric(as.vector(BedFile[,EndPos]))),strand=rep("*",nrow(BedFile)))
      if(ncol(BedFile) > 3){
		elementMetadata(TempRanges_Bed) <- BedFile[,-c(ChrPos,StartPos,EndPos)]
		colnames(elementMetadata(TempRanges_Bed)) <- c("Peak_ID","Score")
	}
      TempRanges_Bed
}

MakeConsensusSet <- function(PeakFileList){

  for(i in 1:length(PeakFileList)){
    if(i == 1){
      ToMerge <- PeakFileList[[i]]
    }else{
      ToMerge <- c(ToMerge,PeakFileList[[i]])

    }
  }
  NonOverlappingSet <- reduce(ToMerge)
  Temp <- do.call(cbind,lapply(PeakFileList,function(x)countOverlaps(NonOverlappingSet,x)))
  Temp[Temp > 1] <- 1
  elementMetadata(NonOverlappingSet) <- Temp
  return(NonOverlappingSet)
}

###################
## Read in peak directory
###################

All_PeakFiles <-  dir(file.path(getwd(),"/Data_For_ChIP_Practical/"),pattern="*_peaks.bed",full.names=T)

PeakFileList <- lapply(All_PeakFiles,Bed2GRanges)
names(PeakFileList) <- gsub("_peaks.bed","",basename(All_PeakFiles))


## Peaks in file 1 which overlap peaks in 2
PeakFileList[[1]][PeakFileList[[1]] %over% PeakFileList[[2]]]
## Peaks in file 2 which overlap peaks in 1
PeakFileList[[2]][PeakFileList[[2]] %over% PeakFileList[[1]]]
## Peaks unique to file 2
PeakFileList[[2]][!PeakFileList[[2]] %over% PeakFileList[[1]]]

###################
## Consensus set
###################
MergedPeakSets <- MakeConsensusSet(PeakFileList)

## Peaks in file 1 which overlap peaks in consensus set
PeakFileList[[1]][PeakFileList[[1]] %over% MergedPeakSets]
## Peaks in file 2 which overlap peaks in consensus set
PeakFileList[[2]][PeakFileList[[2]] %over% MergedPeakSets]
## Peaks unique to file 2
PeakFileList[[2]][!PeakFileList[[2]] %over% MergedPeakSets]


a <- vennCounts(as.data.frame(elementMetadata(MergedPeakSets)))
print(a)

png("Overlap_Diagram.png")
vennDiagram(a)
dev.off()

#####################
## Define peak sets under different conditions
####################

MergedPeakSets <- MergedPeakSets[seqnames(MergedPeakSets) %in% "chr1"]
seqlevels(MergedPeakSets) <- "chr1"

MergedSets_TFonly <- MergedPeakSets[elementMetadata(MergedPeakSets)$TF_1 ==1 & elementMetadata(MergedPeakSets)$TF_2 ==1 & elementMetadata(MergedPeakSets)$CoTF ==0]
MergedSets_TFandCoTF <- MergedPeakSets[elementMetadata(MergedPeakSets)$TF_1 ==1 & elementMetadata(MergedPeakSets)$TF_2 ==1 & elementMetadata(MergedPeakSets)$CoTF ==1]

length(MergedSets_TFonly)
length(MergedSets_TFandCoTF)


#####################
## Use BigWigs to make average signal profile over peaks
####################

bwfileCoTF <- file.path(getwd(),"/Data_For_ChIP_Practical/CoTF.bw")
bwfileTF <- file.path(getwd(),"/Data_For_ChIP_Practical/TF_1.bw")

GetAverageSignalOverRanges <- function(bwfile,selection,Window){
  resizedselection <- resize(selection,fix="center",Window)
  BWSLeft <- BigWigSelection(ranges = resizedselection)
  TempLeft <- import(bwfile,which=resizedselection,selection = BWSLeft)
  CovResLeft <- rep(TempLeft$score,width(TempLeft))
  YouLeft <-   matrix(CovResLeft,ncol=(Window),byrow=T)
  LeftMatrix1 <- apply(YouLeft, MARGIN = 2, FUN = function(X) (X - min(X))/diff(range(X)))
  Leftmeans <- colMeans(YouLeft)
  return(Leftmeans)
}

SignalOverCommon_TF <- GetAverageSignalOverRanges(bwfileTF,MergedSets_TFandCoTF,1000)
SignalOverCommon_CoTF <- GetAverageSignalOverRanges(bwfileCoTF,MergedSets_TFandCoTF,1000)
SignalOverTFOnly_TF <- GetAverageSignalOverRanges(bwfileTF,MergedSets_TFonly,1000)
SignalOverTFOnly_CoTF <- GetAverageSignalOverRanges(bwfileCoTF,MergedSets_TFonly,1000)


png("TrialPlot.png")
par(mfrow=c(1,2))
plot(SignalOverCommon_TF,col="red",ylim=c(0,max(SignalOverTFOnly_CoTF,SignalOverTFOnly_TF,SignalOverCommon_CoTF,SignalOverCommon_TF)*1.2),type="l")
lines(SignalOverCommon_CoTF,col="green")
legend("topright",fill=c("red","green"),legend=c("TF","CoTF"))
plot(SignalOverTFOnly_TF,col="red",ylim=c(0,max(SignalOverTFOnly_CoTF,SignalOverTFOnly_TF,SignalOverCommon_CoTF,SignalOverCommon_TF)*1.2),type="l")
lines(SignalOverTFOnly_CoTF,col="green")
legend("topright",fill=c("red","green"),legend=c("TF","CoTF"))
dev.off()

#############################
## Get counts in Peaks
#############################
SummarisedExperimentList <- summarizeOverlaps(MergedPeakSets,BamsForAnalysis)
assays(SummarisedExperimentList)$counts[1:10,]

elementMetadata(MergedPeakSets) <- cbind(as.data.frame(elementMetadata(MergedPeakSets)),assays(SummarisedExperimentList)$counts)


###########################
### Annotation
###########################


myPeakList1 = as(MergedSets_TFandCoTF,"RangedData") #convert GRanges to RangedData format
data(TSS.human.GRCh37) # Loads gene locations for human genome hg19.
annotatedPeak = annotatePeakInBatch (myPeakList1, AnnotationData = TSS.human.GRCh37) #annotate the peaks
write.table(as.data.frame(annotatedPeak), file = "annotatedPeakList.xls", sep = "\t", row.names = FALSE) #write to spreadsheet
y = annotatedPeak$distancetoFeature[!is.na(annotatedPeak$distancetoFeature) & annotatedPeak$fromOverlappingOrNearest == "NearestStart"] #calculate distance to feature

png("DistanceToNearestTSS.png")
hist(y, xlab = "Distance To Nearest TSS") #plot distance to feature distribution
dev.off()
png("DistributionWindowAroundTSS.png")
hist(y[y<=100 | y>=-100], xlab = "Distance To Nearest TSS") #plot to show the distribution of TF binding around TSS region (-100, +100)
dev.off()

library(org.Hs.eg.db)
enrichedGO = getEnrichedGO (annotatedPeak, orgAnn = "org.Hs.eg.db", maxP =0:01, multiAdj = TRUE, minGOterm = 10, multiAdjMethod = "BH" )
annotatedBDP = peaksNearBDP(myPeakList1, AnnotationData=TSS.human.NCBI37, MaxDistance=5000,PeakLocForDistance = "middle", FeatureLocForDistance = "TSS")


###########################
### Extracting Sequences
###########################


library(BSgenome.Hsapiens.UCSC.hg19)
PeakCentres <- resize(MergedSets_TFandCoTF,fix="center",150)
RangedData_PeakCentres = as(PeakCentres,"RangedData")
Temp <- getAllPeakSequence(RangedData_PeakCentres,upstream=0,downstream=0,genome = Hsapiens) # Extracts corresponding sequences from Human genome
Sequences <- Temp$sequence
names(Sequences) <- paste(Temp$space,start(Temp),end(Temp),sep="_")
SequencesAsXstring <- DNAStringSet(Sequences)
writeXStringSet(SequencesAsXstring,file="test.fasta",width=150)
