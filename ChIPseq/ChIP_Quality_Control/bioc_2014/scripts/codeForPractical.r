################

library(ChIPQC)
load("/data/ChIPQC/ebf1_ChIPQCsample.RData")
QCmetrics(exampleExp)
################

load("/data/ChIPQC/ebf1_FullyLoaded_ChIPQCsample.RData")
QCmetrics(exampleExp)
################

plotCC(exampleExp)
FragmentLengthCrossCoverage(exampleExp)
RelativeCrossCoverage(exampleExp)
################

frip(exampleExp)
plotFrip(exampleExp)
plotFribl(exampleExp)
################

regi(exampleExp)
plotRegi(exampleExp)
################

coveragehistogram(exampleExp)[1:10]
sum(coveragehistogram(exampleExp)[-c(1:19)])
plotCoverageHist(exampleExp)
################

plotSSD(exampleExp)
################


################
################

load("/data/ChIPQC/SampleSheet.RData")
SampleSheet[1:3,]
#resExperiment = ChIPQC(SampleSheet,peaks=peaksFile,annotation="mm9",blacklist=BlackListFile)
load("/data/ChIPQC/BCell_Examples.RData")
QCmetrics(resExperiment)[5:10,]
################

FragmentLengthCrossCoverage(resExperiment)[1:3]
RelativeCrossCoverage(resExperiment)[1:3]
plotCC(resExperiment,colourBy="Tissue",facetBy="Factor",lineBy="Replicate")
################

ccplot <- plotCC(resExperiment,colourBy="Tissue",facetBy="Factor",lineBy="Replicate")
ccplot +facet_wrap(~Factor,scales="free_y")
################

plotFrip(resExperiment)
plotFribl(resExperiment)
################

plotRegi(resExperiment)
regi(resExperiment)["All3utrs",]
plotRegi(resExperiment)+scale_fill_gradient2(low="white",high="red",
         mid="white",midpoint=regi(resExperiment)["All3utrs","Input_Ch12"])
################

plotCoverageHist(resExperiment,facetBy="Factor",colourBy="Tissue",lineBy="Replicate")
################

plotSSD(resExperiment)
################
plotCorHeatmap(resExperiment,facetBy="Factor")

################
################
metrics <- QCmetrics(resExperiment)
metricsMetadata <- data.frame(SampleID=rownames(metrics),RIP =metrics[,"RiP%",drop=T])
plotCC(resExperiment,addMetaData = metricsMetadata,colourBy="RIP")

################

metricsMetadata <- data.frame(SampleID=rownames(metrics),ChIPType=
                                c(rep("Epigenetic Mark",1),
                                  rep("Transcription Factor",1),
                                  rep("Epigenetic Mark",4),
                                  rep("Transcription Factor",4),
                                  rep("Input",2),
                                  rep("Transcription Factor",1),
                                  rep("Epigenetic Mark",2),                    
                                  rep("Transcription Factor",2),                    
                                  rep("Epigenetic Mark",2)))
plotSSD(resExperiment,addMetaData = metricsMetadata,facetBy="ChIPType")
################

DNAsePeaks <- peaks(QCsample(resExperiment,"DNAse"))
H3K4me3Peaks <- peaks(QCsample(resExperiment,"H3K4me3_IkPos"))
H3K9acPeaks <- peaks(QCsample(resExperiment,"H3K9ac_IkPos"))
RNAPol2Peaks <- peaks(QCsample(resExperiment,"RNAPol2"))

customAnnotation <- list(version="custom",
                         DNAse=DNAsePeaks,
                         RNAPol2=RNAPol2Peaks,
                         H3K4me3=H3K4me3Peaks,
                         H3K9ac=H3K9acPeaks)
#bamFile <- "/data/ChIPQC/Chr11_Ebf1DupMarked.bam
bamFile <- "/data/ChIPQC/Chr11_Ebf1DupMarked.bam"
#exampleExpCA = ChIPQCsample(bamFile,peaks=NULL,annotation=customAnnotation,chromosomes="chr11")
load("/data/ChIPQC/exampleCustomAnnotation.RData")
plotRegi(exampleExpCA)
################

ccplot <- plotCC(resExperiment)
ccplot$layers <- ccplot$layers[1]
ccplot
################

plotFribl(resExperiment)+ylim(0,15)
################

plotCoverageHist(resExperiment)+xlim(0,250)

################
################
################
data(tamoxifen_QC)
QCmetrics(tamoxifen)
## Some things to try
plotCC(tamoxifen)
plotCC(tamoxifen,colourBy="Tissue")+facet_wrap(~Factor,scales="free")
plotFrip(tamoxifen)
plotFribl(tamoxifen)+ylim(0,20)
plotCoverageHist(tamoxifen)+xlim(0,1000)
plotSSD(tamoxifen)+xlim(0,10)

