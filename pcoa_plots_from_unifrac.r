#!/usr/bin/env Rscript
options(error=recover)

#this script prints out PDF pcoa plots and distance matrices, given an OTU table, phylogenetic tree, and metadata

source("UniFrac.r")
library(ape)
library(phangorn)
library(vegan)

otu.tab <- read.table("data/td_OTU_tag_mapped.txt", header=T, sep="\t", row.names=1, comment.char="", check.names=FALSE)

otu.tab <- t(as.matrix(otu.tab))

#sort taxa from most to least abundant
taxaOrder <- rev(order(apply(otu.tab,2,sum)))
otu.tab <- otu.tab[,taxaOrder]

# read and root tree (rooted tree is required)
tree <- read.tree("data/OTU_seqs.tre")
tree <- midpoint(tree)

# read metadata
MyMeta<- read.table("data/metadata.txt", header=T, sep="\t", row.names=1, comment.char="", check.names=FALSE)

# clean up metadata
metadata <- MyMeta
# first column has nothing except an N in the Total row
metadata <- metadata[,c(2:ncol(metadata))]
# remove Total row
metadata <- metadata[c(1:(nrow(metadata)-1)),]
# make Sex be factor of either M or F (there's a 316 in there from Total)
metadata$Sex <- as.character(metadata$Sex)
metadata$Sex <- as.factor(metadata$Sex)
# make smoker be factor of N Q or Y (there's a 316 in there from Total)
metadata$smoker <- as.character(metadata$smoker)
metadata$smoker <- as.factor(metadata$smoker)

# filter OTU table and metadata so that only samples which appear in both are retained
otu_indicies <- match(rownames(metadata),rownames(otu.tab))
# this will remove [1] "mix_167_20", a contaminated sample
otu_indicies <- otu_indicies[!is.na(otu_indicies)]
otu.tab <- otu.tab[otu_indicies,]
# put metadata samples in the same order as otu.tab
metadata <- metadata[match(rownames(otu.tab),rownames(metadata)),]

# 863 OTUs in total - remove anything less than 1% abundant in all samples
otu.tab.sum <- apply(otu.tab,1,sum)
otu.tab.sum.threshhold <- otu.tab.sum * 0.01
otu.tab.filter <- apply(otu.tab,2, function(x) { return(length(which(x > otu.tab.sum.threshhold))) } )
otu.tab <- otu.tab[,which(otu.tab.filter > 0)]

# remove low abundance OTUs from tree
absent <- tree$tip.label[!(tree$tip.label %in% colnames(otu.tab))]
if (length(absent) != 0) {
		tree <- drop.tip(tree, absent)
}

#rarefy data for unweighted unifrac to min count per sample (9654)
otu.tab.rarefy <- rrarefy(otu.tab, min(apply(otu.tab,1,sum)))

#calculate distance matrix
unweighted <- getDistanceMatrix(otu.tab.rarefy,tree,method="unweighted",verbose=TRUE)
all.dist.mat <- getDistanceMatrix(otu.tab,tree,method="all",verbose=TRUE)

weighted <- all.dist.mat[["weighted"]]
information <- all.dist.mat[["information"]]
ratio_no_log <- all.dist.mat[["ratio_no_log"]]

#output distance matrices
write.table(unweighted,file="unweighted_distance_matrix.txt",sep="\t",quote=FALSE)
write.table(weighted,file="weighted_distance_matrix.txt",sep="\t",quote=FALSE)
write.table(information,file="information_distance_matrix.txt",sep="\t",quote=FALSE)
write.table(ratio_no_log,file="ratio_normalize_distance_matrix.txt",sep="\t",quote=FALSE)

colnames(metadata) <- gsub(" ",".",colnames(metadata))
residuals <- metadata$Standardized.Residual
groups <- rep("NA",length(residuals))
groups[which(residuals < -1)] <- "Protected"
groups[which(residuals >= -1 & residuals < 1)] <- "Explained"
groups[which(residuals >= 1)] <- "Unexplained"
groups <- as.factor(groups)

unweighted.pcoa <- pcoa(unweighted)
weighted.pcoa <- pcoa(weighted)
information.pcoa <- pcoa(information)
ratio_no_log.pcoa <- pcoa(ratio_no_log)

#function to get variance explained for the PCOA component labels
getVarExplained <- function(vector) {
	rawVarEx <- apply(vector,2,function(x) sd(x)*sd(x))
	totalVarExplained <- sum(rawVarEx)
	varEx <- rawVarEx/totalVarExplained
	return(varEx)
}

unweighted.varEx <- getVarExplained(unweighted.pcoa$vectors)
weighted.varEx <- getVarExplained(weighted.pcoa$vectors)
information.varEx <- getVarExplained(information.pcoa$vectors)
ratio_no_log.varEx <- getVarExplained(ratio_no_log.pcoa$vectors)

#get vector version of distance matrices for correlation plots below
unweighted.vector <- unlist(unweighted[lower.tri(unweighted,diag=TRUE)])
weighted.vector <- unlist(weighted[lower.tri(weighted,diag=TRUE)])
information.vector <- unlist(information[lower.tri(information,diag=TRUE)])
ratio_no_log.vector <- unlist(ratio_no_log[lower.tri(ratio_no_log,diag=TRUE)])

plot.as.factor <- function(unweighted.pcoa, weighted.pcoa, information.pcoa, ratio.pcoa, factor, name, my.palette) {
	unweighted.varEx <- getVarExplained(unweighted.pcoa$vectors)
	weighted.varEx <- getVarExplained(weighted.pcoa$vectors)
	information.varEx <- getVarExplained(information.pcoa$vectors)
	ratio_no_log.varEx <- getVarExplained(ratio_no_log.pcoa$vectors)
	
	factor <- as.factor(factor)
	palette(my.palette)
	par(mar=c(5.1, 5.1, 5.1, 7.1),xpd=TRUE)
	plot(unweighted.pcoa$vectors[,1],unweighted.pcoa$vectors[,2], col=factor,main=paste("Unweighted UniFrac",name,sep="\n"),xlab=paste("First Component", round(unweighted.varEx[1],digits=3),"variance explained"),ylab=paste("Second Component", round(unweighted.varEx[2],digits=3),"variance explained"),pch=19,cex.lab=1.4,cex.main=2)
	legend("topright", levels(factor), pch=rep(19,length(factor)), col=palette(), xpd=NA, inset=c(-0.3,0))
	plot(weighted.pcoa$vectors[,1],weighted.pcoa$vectors[,2], col=factor,main=paste("Weighted UniFrac",name,sep="\n"),xlab=paste("First Component", round(weighted.varEx[1],digits=3),"variance explained"),ylab=paste("Second Component", round(weighted.varEx[2],digits=3),"variance explained"),pch=19,cex.lab=1.4,cex.main=2)
	legend("topright", levels(factor), pch=rep(19,length(factor)), col=palette(), xpd=NA, inset=c(-0.3,0))
	plot(information.pcoa$vectors[,1],information.pcoa$vectors[,2], col=factor,main=paste("Information UniFrac",name,sep="\n"),xlab=paste("First Component", round(information.varEx[1],digits=3),"variance explained"),ylab=paste("Second Component", round(information.varEx[2],digits=3),"variance explained"),pch=19,cex.lab=1.4,cex.main=2)
	legend("topright", levels(factor), pch=rep(19,length(factor)), col=palette(), xpd=NA, inset=c(-0.3,0))
	plot(ratio_no_log.pcoa$vectors[,1],ratio_no_log.pcoa$vectors[,2], col=factor,main=paste("Ratio UniFrac",name,sep="\n"),xlab=paste("First Component", round(ratio_no_log.varEx[1],digits=3),"variance explained"),ylab=paste("Second Component", round(ratio_no_log.varEx[2],digits=3),"variance explained"),pch=19,cex.lab=1.4,cex.main=2)	
	legend("topright", levels(factor), pch=rep(19,length(factor)), col=palette(), xpd=NA, inset=c(-0.3,0))
}

plot.as.number <- function(unweighted.pcoa, weighted.pcoa, information.pcoa, ratio.pcoa, numbers, name, my.palette) {
	unweighted.varEx <- getVarExplained(unweighted.pcoa$vectors)
	weighted.varEx <- getVarExplained(weighted.pcoa$vectors)
	information.varEx <- getVarExplained(information.pcoa$vectors)
	ratio_no_log.varEx <- getVarExplained(ratio_no_log.pcoa$vectors)

	summary(numbers)
	groups <- rep("NA",length(numbers))
	groups.quart.1 <- summary(numbers)[2]
	groups.quart.3 <- summary(numbers)[5]
	groups[which(numbers <= groups.quart.1)] <- "1st Quartile"
	groups[which((numbers > groups.quart.1) & (numbers < groups.quart.3))] <- "Middle Half"
	groups[which(numbers >= groups.quart.3)] <- "3rd Quartile"
	groups <- as.factor(groups)
	palette(my.palette)
	par(mar=c(5.1, 5.1, 5.1, 7.1),xpd=TRUE)
	plot(unweighted.pcoa$vectors[,1],unweighted.pcoa$vectors[,2], col=groups,main=paste("Unweighted UniFrac",name,sep="\n"),xlab=paste("First Component", round(unweighted.varEx[1],digits=3),"variance explained"),ylab=paste("Second Component", round(unweighted.varEx[2],digits=3),"variance explained"),pch=19,cex.lab=1.4,cex.main=2)
	legend("topright", levels(groups), pch=rep(19,length(groups)), col=palette(), xpd=NA, inset=c(-0.3,0))
	plot(weighted.pcoa$vectors[,1],weighted.pcoa$vectors[,2], col=groups,main=paste("Weighted UniFrac",name,sep="\n"),xlab=paste("First Component", round(weighted.varEx[1],digits=3),"variance explained"),ylab=paste("Second Component", round(weighted.varEx[2],digits=3),"variance explained"),pch=19,cex.lab=1.4,cex.main=2)
	legend("topright", levels(groups), pch=rep(19,length(groups)), col=palette(), xpd=NA, inset=c(-0.3,0))
	plot(information.pcoa$vectors[,1],information.pcoa$vectors[,2], col=groups,main=paste("Information UniFrac",name,sep="\n"),xlab=paste("First Component", round(information.varEx[1],digits=3),"variance explained"),ylab=paste("Second Component", round(information.varEx[2],digits=3),"variance explained"),pch=19,cex.lab=1.4,cex.main=2)
	legend("topright", levels(groups), pch=rep(19,length(groups)), col=palette(), xpd=NA, inset=c(-0.3,0))
	plot(ratio_no_log.pcoa$vectors[,1],ratio_no_log.pcoa$vectors[,2], col=groups,main=paste("Ratio UniFrac",name,sep="\n"),xlab=paste("First Component", round(ratio_no_log.varEx[1],digits=3),"variance explained"),ylab=paste("Second Component", round(ratio_no_log.varEx[2],digits=3),"variance explained"),pch=19,cex.lab=1.4,cex.main=2)
	legend("topright", levels(groups), pch=rep(19,length(groups)), col=palette(), xpd=NA, inset=c(-0.3,0))
}

original.palette <- palette()
pdf("pcoa_plots.pdf")

plot.as.factor(unweighted.pcoa, weighted.pcoa, information.pcoa, ratio_no_log.pcoa, groups, "principal coordinate analysis", c("purple","blue","red"))

plot.as.number(unweighted.pcoa, weighted.pcoa, information.pcoa, ratio_no_log.pcoa, metadata$TPAmm2, "TPA mm2", c("blue","red","purple"))

plot.as.number(unweighted.pcoa, weighted.pcoa, information.pcoa, ratio_no_log.pcoa, metadata$Age, "Age", c("blue","red","purple"))

plot.as.factor(unweighted.pcoa, weighted.pcoa, information.pcoa, ratio_no_log.pcoa, metadata$Sex, "Sex", c("blue","red","purple"))

plot.as.factor(unweighted.pcoa, weighted.pcoa, information.pcoa, ratio_no_log.pcoa, metadata$smoker, "Smoking", c("blue","purple","red"))

plot.as.number(unweighted.pcoa, weighted.pcoa, information.pcoa, ratio_no_log.pcoa, metadata$PackYears, "Pack Years", c("blue","red","purple"))

plot.as.number(unweighted.pcoa, weighted.pcoa, information.pcoa, ratio_no_log.pcoa, metadata$BPSys, "Systolic Blood Pressure", c("blue","red","purple"))

plot.as.number(unweighted.pcoa, weighted.pcoa, information.pcoa, ratio_no_log.pcoa, metadata$BPDias, "Diastolic Blood Pressure", c("blue","red","purple"))

plot.as.number(unweighted.pcoa, weighted.pcoa, information.pcoa, ratio_no_log.pcoa, metadata$chol, "Cholesterol", c("blue","red","purple"))

plot.as.number(unweighted.pcoa, weighted.pcoa, information.pcoa, ratio_no_log.pcoa, metadata$trig, "Triglycerides", c("blue","red","purple"))

plot.as.number(unweighted.pcoa, weighted.pcoa, information.pcoa, ratio_no_log.pcoa, metadata$hdl, "HDL", c("blue","red","purple"))

plot.as.number(unweighted.pcoa, weighted.pcoa, information.pcoa, ratio_no_log.pcoa, metadata$ldl, "LDL", c("blue","red","purple"))

dev.off()

residuals <- metadata$Standardized.Residual
top <- which(residuals >= 2)
bottom <- which(residuals <= -2)
otu.tab.rarefy.extreme <- otu.tab.rarefy[c(top,bottom),]
otu.tab.extreme <- otu.tab[c(top,bottom),]
groups <- c(rep("Unexplained",length(top)),rep("Protected",length(bottom)))
groups <- as.factor(groups)
# all OTUs have at least 1 count :)

#calculate distance matrix
unweighted <- getDistanceMatrix(otu.tab.rarefy.extreme,tree,method="unweighted",verbose=TRUE)
all.dist.mat <- getDistanceMatrix(otu.tab.extreme,tree,method="all",verbose=TRUE)

weighted <- all.dist.mat[["weighted"]]
information <- all.dist.mat[["information"]]
ratio_no_log <- all.dist.mat[["ratio_no_log"]]

#output distance matrices
write.table(unweighted,file="unweighted_distance_matrix_extreme_deciles.txt",sep="\t",quote=FALSE)
write.table(weighted,file="weighted_distance_matrix_extreme_deciles.txt",sep="\t",quote=FALSE)
write.table(information,file="information_distance_matrix_extreme_deciles.txt",sep="\t",quote=FALSE)
write.table(ratio_no_log,file="ratio_normalize_distance_matrix_extreme_deciles.txt",sep="\t",quote=FALSE)

unweighted.pcoa <- pcoa(unweighted)
weighted.pcoa <- pcoa(weighted)
information.pcoa <- pcoa(information)
ratio_no_log.pcoa <- pcoa(ratio_no_log)

pdf("pcoa_plots_extreme_residuals.pdf")

plot.as.factor(unweighted.pcoa, weighted.pcoa, information.pcoa, ratio_no_log.pcoa, groups, "principal coordinate analysis", c("purple","blue","red"))

plot.as.number(unweighted.pcoa, weighted.pcoa, information.pcoa, ratio_no_log.pcoa, metadata$TPAmm2, "TPA mm2", c("blue","red","purple"))

plot.as.number(unweighted.pcoa, weighted.pcoa, information.pcoa, ratio_no_log.pcoa, metadata$Age, "Age", c("blue","red","purple"))

plot.as.factor(unweighted.pcoa, weighted.pcoa, information.pcoa, ratio_no_log.pcoa, metadata$Sex, "Sex", c("blue","red","purple"))

plot.as.factor(unweighted.pcoa, weighted.pcoa, information.pcoa, ratio_no_log.pcoa, metadata$smoker, "Smoking", c("blue","purple","red"))

plot.as.number(unweighted.pcoa, weighted.pcoa, information.pcoa, ratio_no_log.pcoa, metadata$PackYears, "Pack Years", c("blue","red","purple"))

plot.as.number(unweighted.pcoa, weighted.pcoa, information.pcoa, ratio_no_log.pcoa, metadata$BPSys, "Systolic Blood Pressure", c("blue","red","purple"))

plot.as.number(unweighted.pcoa, weighted.pcoa, information.pcoa, ratio_no_log.pcoa, metadata$BPDias, "Diastolic Blood Pressure", c("blue","red","purple"))

plot.as.number(unweighted.pcoa, weighted.pcoa, information.pcoa, ratio_no_log.pcoa, metadata$chol, "Cholesterol", c("blue","red","purple"))

plot.as.number(unweighted.pcoa, weighted.pcoa, information.pcoa, ratio_no_log.pcoa, metadata$trig, "Triglycerides", c("blue","red","purple"))

plot.as.number(unweighted.pcoa, weighted.pcoa, information.pcoa, ratio_no_log.pcoa, metadata$hdl, "HDL", c("blue","red","purple"))

plot.as.number(unweighted.pcoa, weighted.pcoa, information.pcoa, ratio_no_log.pcoa, metadata$ldl, "LDL", c("blue","red","purple"))

dev.off()
