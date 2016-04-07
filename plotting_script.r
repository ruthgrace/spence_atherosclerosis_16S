library(zCompositions)
library(randomcoloR)
library(compositions)
library(ALDEx2)
library(stringr)

metagenomic_samples <- c("CL_119", "CL_139", "CL_141", "CL_144", "CL_160", "CL_165", "CL_166", "CL_169", "CL_173", "CL_177", "HLD_100", "HLD_102", "HLD_111", "HLD_112", "HLD_23", "HLD_28", "HLD_47", "HLD_72", "HLD_80", "HLD_85")

# read metadata for 16S samples
MyMeta<- read.table("../exponentUnifrac/data/nash_data/metadata.txt", header=T, sep="\t", row.names=1, comment.char="", check.names=FALSE)
metadata <- MyMeta[grepl("a$",rownames(MyMeta)),]
rownames(metadata) <- gsub("a$","",rownames(metadata))
samples <- str_extract(rownames(metadata), "^[A-Z]*-[0-9]*")
samples <- gsub("-","_",samples)
unique.samples <- unique(samples)
metadata <- metadata[match(unique.samples,samples),]
rownames(metadata) <- unique.samples

# get metaphlan data
d <- read.table("data/summary_all_count.txt", header=T, row.names=1, sep="\t",quote="",comment.char="",stringsAsFactors=FALSE)

d <- d[,(grepl("CL*", colnames(d)) | grepl("HLD*", colnames(d)))]

# for some reason the unclassified species has a bunch of dashes instead of numbers
d <- d[which(rownames(d)!="s__unclassified"),]
d.rownames <- rownames(d)
d <- apply(d,2,function(x) as.numeric(x))
rownames(d) <- d.rownames

# remove all features with zero counts for all samples
d.sum <- apply(d,1,sum)

d <- d[which(d.sum>0),]

original.data <- d

rownames(d) <- gsub("s__","",rownames(d))
rownames(d) <- sub("_"," ",rownames(d))

# get genus level for metaphlan
species <- rownames(d)
genus <- str_extract(species,"^[A-Za-z]*")

d.genus <- aggregate(d,list(genus),sum)
rownames(d.genus) <- d.genus$Group.1
d.genus <- d.genus[,c(2:ncol(d.genus))]
colnames(d.genus) <- str_extract(colnames(d.genus),"^[A-Z]*_[0-9]*")

# get genus level for 16S
otu.tab <- read.table("data/summed_data_baseline_only_with_taxonomy.txt", header=T, sep="\t", row.names=1, comment.char="", check.names=FALSE)

colnames(otu.tab)[c(1:(ncol(otu.tab)-1))] <- gsub("-","_",colnames(otu.tab)[c(1:(ncol(otu.tab)-1))])

colnames(otu.tab)[c(1:(ncol(otu.tab)-1))] <- str_extract(colnames(otu.tab)[c(1:(ncol(otu.tab)-1))],"[A-Z]+_[0-9]+")

taxonomy <- otu.tab$taxonomy

otu.tab <- otu.tab[,c(1:(ncol(otu.tab)-1))]
otu.genus <- c(as.character(taxonomy))

for (i in c(1:length(taxonomy))) {
  otu.genus[i] <- strsplit(otu.genus[i],c(";"))[[1]][6]
}

otu.tab.genus <- aggregate(otu.tab,list(otu.genus),sum)
rownames(otu.tab.genus) <- otu.tab.genus$Group.1
otu.tab.genus <- otu.tab.genus[,c(2:ncol(otu.tab.genus))]

remove <- c("Incertae_Sedis", "unclassified")

print(paste(sum(otu.tab.genus["Incertae_Sedis",]), " Incertae_Sedis out of",sum(otu.tab.genus),"total counts"))
# [1] "1156172  Incertae_Sedis out of 7472468 total counts"
# > 1156172/7472468
# [1] 0.1547242

print(paste(sum(otu.tab.genus["unclassified",]), " unclassified out of",sum(otu.tab.genus),"total counts"))
# [1] "937127  unclassified out of 7472468 total counts"
# > 937127/7472468
# [1] 0.1254106

otu.tab.genus <- otu.tab.genus[which(! rownames(otu.tab.genus) %in% remove),]

# conditions: Originally 0 meant steatohepatosis, and 1 meant NASH
groups <- metadata$SSvsNASH[match(colnames(otu.tab.genus),rownames(metadata))]
originalgroups <- groups

# Make healthy represented by 0, SS by 1, NASH by 2
groups <- groups + 1;
groups[which(is.na(groups))] <- 0

# make healthy 1, ss 2, nash 3 (healthy metagenomic will be 0 and nash metagenomic will be 4)
groups <- groups + 1

# mark healthy samples selected for metagenomic study
groups[which(colnames(otu.tab.genus) %in% metagenomic_samples & groups == 1)] <- 0

# mark nash samples selected for metagenomic study
groups[which(colnames(otu.tab.genus) %in% metagenomic_samples & groups == 3)] <- 4

groups[which(groups == 0)] <- "Healthy Metagenomic"
groups[which(groups == 1)] <- "Healthy"
groups[which(groups == 2)] <- "SS"
groups[which(groups == 3)] <- "NASH"
groups[which(groups == 4)] <- "NASH Metagenomic"

groups <- as.factor(groups)


sample.sum <- apply(d,2,sum)
one.percent <- sample.sum*0.01

# adjust zeros
d.adj.zero <- t(cmultRepl(t(d),method="CZM"))
d.genus.adj.zero <- t(cmultRepl(t(d.genus),method="CZM"))
otu.tab.genus.adj.zero <- t(cmultRepl(t(otu.tab.genus),method="CZM"))

filter <- apply(d,1,function(x) length(which(x > one.percent)))
d.filter <- d.adj.zero[which(filter > 0),]

d.filter.counts <- d[which(filter>0),]

d.adj.zero <- d.adj.zero[order(apply(d.adj.zero,1,sum),decreasing=TRUE),]
d.filter <- d.filter[order(apply(d.filter,1,sum),decreasing=TRUE),]
d.genus.adj.zero <- d.genus.adj.zero[order(apply(d.genus.adj.zero,1,sum),decreasing=TRUE),]
otu.tab.genus.adj.zero <- otu.tab.genus.adj.zero[order(apply(otu.tab.genus.adj.zero,1,sum),decreasing=TRUE),]

d.names <- rownames(d.adj.zero)
d.filter.names <- rownames(d.filter)
d.genus.names <- rownames(d.genus.adj.zero)
otu.tab.genus.names <- rownames(otu.tab.genus.adj.zero)

taxa.col <- data.frame(as.character(rownames(d)),rownames(d))
colnames(taxa.col) <- c("taxon","color")
taxa.col[,2] <- distinctColorPalette(length(taxa.col[,2]))

taxa.filter.col <- data.frame(as.character(rownames(d.filter)),rownames(d.filter))
colnames(taxa.filter.col) <- c("taxon","color")
taxa.filter.col[,2] <- taxa.col[match(taxa.filter.col[,1],taxa.col[,1]),2]

all.genus <- unique(c(rownames(d.genus),rownames(otu.tab.genus)))
all.genus.colors <- distinctColorPalette(length(all.genus))

taxa.d.genus.col <- data.frame(rownames(d.genus),rownames(d.genus))
colnames(taxa.d.genus.col) <- c("taxon","color")
taxa.d.genus.col[,2] <- all.genus.colors[match(rownames(d.genus),all.genus)]

taxa.otu.tab.genus.col <- data.frame(rownames(otu.tab.genus),rownames(otu.tab.genus))
colnames(taxa.otu.tab.genus.col) <- c("taxon","color")
taxa.otu.tab.genus.col[,2] <- all.genus.colors[match(rownames(otu.tab.genus),all.genus)]


d.prop <- apply(d.adj.zero,2,function(x){x/sum(x)})
d.filter.prop <- apply(d.filter,2,function(x) {x/sum(x)})
d.genus.prop <- apply(d.genus.adj.zero, 2,function(x) {x/sum(x)})
otu.tab.genus.prop <- apply(otu.tab.genus.adj.zero, 2,function(x) {x/sum(x)})


d.clr <- t(apply(d.prop,2,function(x){log(x) - mean(log(x))}))
d.filter.clr <- t(apply(d.filter.prop,2,function(x){log(x) - mean(log(x))}))
d.genus.clr <- t(apply(d.genus.prop,2,function(x){log(x) - mean(log(x))}))
otu.tab.genus.clr <- t(apply(otu.tab.genus.prop,2,function(x){log(x) - mean(log(x))}))


d.pcx <- prcomp(d.clr)
d.filter.pcx <- prcomp(d.filter.clr)
d.genus.pcx <- prcomp(d.genus.clr)
otu.tab.genus.pcx <- prcomp(otu.tab.genus.clr)

conds <- data.frame(c(rep("NASH",10),rep("Healthy",10)))
colnames(conds) <- "cond"

palette=palette(c(rgb(1,0,0,0.6), rgb(0,0,1,0.6), rgb(0,1,1,0.6)))

pdf("biplots.pdf")

layout(matrix(c(1,2),1,2, byrow=T), widths=c(6,2), heights=c(8,3))
par(mgp=c(2,0.5,0))
# make a covariance biplot of the data with compositions function
coloredBiplot(d.pcx, cex=c(0.6, 0.6),
arrow.len=0.05,
xlab=paste("PC1 ", round (sum(d.pcx$sdev[1]^2)/mvar(d.clr),3), sep=""),
ylab=paste("PC2 ", round (sum(d.pcx$sdev[2]^2)/mvar(d.clr),3), sep=""),
xlabs.col=c(rep("red",10),rep("black",10)),
expand=0.8,var.axes=FALSE, scale=1, main="Biplot")
barplot(d.pcx$sdev^2/mvar(d.clr),  ylab="variance explained", xlab="Component", main="Scree plot") # scree plot

coloredBiplot(d.filter.pcx, cex=c(0.6, 0.6),
arrow.len=0.05,
xlab=paste("PC1 ", round (sum(d.filter.pcx$sdev[1]^2)/mvar(d.filter.clr),3), sep=""),
ylab=paste("PC2 ", round (sum(d.filter.pcx$sdev[2]^2)/mvar(d.filter.clr),3), sep=""),
xlabs.col=c(rep("red",10),rep("black",10)),
expand=0.8,var.axes=FALSE, scale=1, main="Biplot")
barplot(d.filter.pcx$sdev^2/mvar(d.filter.clr),  ylab="variance explained", xlab="Component", main="Scree plot") # scree plot

coloredBiplot(d.genus.pcx, cex=c(0.6, 0.6),
arrow.len=0.05,
xlab=paste("PC1 ", round (sum(d.genus.pcx$sdev[1]^2)/mvar(d.genus.clr),3), sep=""),
ylab=paste("PC2 ", round (sum(d.genus.pcx$sdev[2]^2)/mvar(d.genus.clr),3), sep=""),
xlabs.col=c(rep("red",10),rep("black",10)),
expand=0.8,var.axes=FALSE, scale=1, main="Biplot")
barplot(d.genus.pcx$sdev^2/mvar(d.genus.clr),  ylab="variance explained", xlab="Component", main="Scree plot") # scree plot

coloredBiplot(otu.tab.genus.pcx, cex=c(0.6, 0.6),
arrow.len=0.05,
xlab=paste("PC1 ", round (sum(otu.tab.genus.pcx$sdev[1]^2)/mvar(otu.tab.genus.clr),3), sep=""),
ylab=paste("PC2 ", round (sum(otu.tab.genus.pcx$sdev[2]^2)/mvar(otu.tab.genus.clr),3), sep=""),
xlabs.col=c(rep("red",10),rep("black",10)),
expand=0.8,var.axes=FALSE, scale=1, main="Biplot")
barplot(otu.tab.genus.pcx$sdev^2/mvar(otu.tab.genus.clr),  ylab="variance explained", xlab="Component", main="Scree plot") # scree plot

dev.off()


# generate the distance matrix
d.dist <- dist(d.clr, method="euclidian")
d.filter.dist <- dist(d.filter.clr, method="euclidian")
d.genus.dist <- dist(d.genus.clr, method="euclidian")
otu.tab.genus.dist <- dist(otu.tab.genus.clr, method="euclidian")

# add condition onto labels of 16S hclust data
otu.tab.conditions <- as.character(groups)
otu.tab.conditions <- gsub(" ","_",otu.tab.conditions)
attributes(otu.tab.genus.dist)$Labels <- paste(otu.tab.conditions, attributes(otu.tab.genus.dist)$Labels, sep="_")

# cluster the data
d.hc <- hclust(d.dist, method="ward.D2")
d.filter.hc <- hclust(d.filter.dist, method="ward.D2")
d.genus.hc <- hclust(d.genus.dist, method="ward.D2")
otu.tab.genus.hc <- hclust(otu.tab.genus.dist, method="ward.D2")

# now re-order the data to plot the barplot in the same order
d.order <- d.adj.zero[,d.hc$order]
d.filter.order <- d.filter[,d.filter.hc$order]
d.genus.order <- d.genus.adj.zero[,d.genus.hc$order]
otu.tab.genus.order <- otu.tab.genus.adj.zero[,otu.tab.genus.hc$order]

d.acomp <- acomp(t(d.order))
d.filter.acomp <- acomp(t(d.filter.order))
d.genus.acomp <- acomp(t(d.genus.order))
otu.tab.genus.acomp <- acomp(t(otu.tab.genus.order))

pdf("dendogram_barplot.pdf")

layout(matrix(c(1,3,2,3),2,2, byrow=T), widths=c(8,10), height=c(4,4))
par(mar=c(2,1,1,1)+0.1)
# plot the dendrogram
plot(d.hc, cex=0.6)
# plot the barplot below
barplot(d.acomp, legend.text=F, col=as.character(taxa.col[,2]), axisnames=F, border=NA, xpd=T)
par(mar=c(0,1,1,1)+0.1)
# and the legend
plot(1,2, pch = 1, lty = 1, ylim=c(-20,20), type = "n", axes = FALSE, ann = FALSE)
legend(x="center", legend=d.names, col=as.character(taxa.col[,2]), lwd=5, cex=.3, border=NULL,ncol=3)

layout(matrix(c(1,3,2,3),2,2, byrow=T), widths=c(6,4), height=c(4,4))
# plot the dendrogram
plot(d.filter.hc, cex=0.6)
# plot the barplot below
barplot(d.filter.acomp, legend.text=F, col=as.character(taxa.filter.col[,2]), axisnames=F, border=NA, xpd=T)
par(mar=c(0,1,1,1)+0.1)
# and the legend
plot(1,2, pch = 1, lty = 1, ylim=c(-20,20), type = "n", axes = FALSE, ann = FALSE)
legend(x="center", legend=d.filter.names, col=as.character(taxa.filter.col[,2]), lwd=5, cex=.5, border=NULL)

layout(matrix(c(1,3,2,3),2,2, byrow=T), widths=c(8,6), height=c(4,4))
# plot the dendrogram
plot(d.genus.hc, cex=0.6)
# plot the barplot below
barplot(d.genus.acomp, legend.text=F, col=as.character(taxa.d.genus.col[,2]), axisnames=F, border=NA, xpd=T)
par(mar=c(0,1,1,1)+0.1)
# and the legend
plot(1,2, pch = 1, lty = 1, ylim=c(-20,20), type = "n", axes = FALSE, ann = FALSE)
legend(x="center", legend=d.genus.names, col=as.character(taxa.d.genus.col[,2]), lwd=5, cex=.5, border=NULL,ncol=2)

layout(matrix(c(1,3,2,3),2,2, byrow=T), widths=c(10,6), height=c(4,4))
# plot the dendrogram
plot(otu.tab.genus.hc, cex=0.4, hang=-1)
# plot the barplot below
barplot(otu.tab.genus.acomp, legend.text=F, col=as.character(taxa.otu.tab.genus.col[,2]), axisnames=F, border=NA, xpd=T)
par(mar=c(0,1,1,1)+0.1)
# and the legend
plot(1,2, pch = 1, lty = 1, ylim=c(-20,20), type = "n", axes = FALSE, ann = FALSE)
legend(x="center", legend=otu.tab.genus.names, col=as.character(taxa.otu.tab.genus.col[,2]), lwd=5, cex=.5, border=NULL,ncol=2)

dev.off()

# generate the dataset by making a data frame of
d.h <- colnames(d)[grep("HLD*", colnames(d))] # Before samples
d.n <- colnames(d)[grep("CL*", colnames(d))] # After samples
d.aldex <- data.frame(d[,d.h], d[,d.n]) # make a data frame
# make the vector of set membership in the same order as
conds.aldex <- c(rep("Healthy", 10), rep("NASH", 10))
# generate 128 Dirichlet Monte-Carlo replicates
x <- aldex.clr(d.aldex, mc.samples=128, verbose=FALSE)
## [1] "operating in serial mode"
# calculate p values for each replicate and report the mean
x.t <- aldex.ttest(x, conds.aldex)
# calculate mean effect sizes
x.e <- aldex.effect(x, conds.aldex, verbose=FALSE)
## [1] "operating in serial mode"
# save it all in a data frame
x.all <- data.frame(x.e,x.t)

# generate the dataset by making a data frame of
d.filter.h <- colnames(d.filter.counts)[grep("HLD*", colnames(d.filter.counts))] # Before samples
d.filter.n <- colnames(d.filter.counts)[grep("CL*", colnames(d.filter.counts))] # After samples
d.filter.aldex <- data.frame(d.filter.counts[,d.filter.h], d.filter.counts[,d.filter.n]) # make a data frame
# make the vector of set membership in the same order as
conds.aldex <- c(rep("Healthy", 10), rep("NASH", 10))
# generate 128 Dirichlet Monte-Carlo replicates
x.filter <- aldex.clr(d.filter.aldex, mc.samples=128, verbose=FALSE)
## [1] "operating in serial mode"
# calculate p values for each replicate and report the mean
x.filter.t <- aldex.ttest(x.filter, conds.aldex)
# calculate mean effect sizes
x.filter.e <- aldex.effect(x.filter, conds.aldex, verbose=FALSE)
## [1] "operating in serial mode"
# save it all in a data frame
x.filter.all <- data.frame(x.filter.e,x.filter.t)

# generate the dataset by making a data frame of
d.genus.h <- colnames(d.genus)[grep("HLD*", colnames(d.genus))] # Before samples
d.genus.n <- colnames(d.genus)[grep("CL*", colnames(d.genus))] # After samples
d.genus.aldex <- data.frame(d.genus[,d.genus.h], d.genus[,d.genus.n]) # make a data frame
# make the vector of set membership in the same order as
conds.aldex <- c(rep("Healthy", 10), rep("NASH", 10))
# generate 128 Dirichlet Monte-Carlo replicates
x.genus <- aldex.clr(d.genus.aldex, mc.samples=128, verbose=FALSE)
## [1] "operating in serial mode"
# calculate p values for each replicate and report the mean
x.genus.t <- aldex.ttest(x.genus, conds.aldex)
# calculate mean effect sizes
x.genus.e <- aldex.effect(x.genus, conds.aldex, verbose=FALSE)
## [1] "operating in serial mode"
# save it all in a data frame
x.genus.all <- data.frame(x.genus.e,x.genus.t)

# generate the dataset by making a data frame of
otu.tab.genus.h <- colnames(otu.tab.genus)[grepl("^Healthy*",as.character(groups))] # Before samples
otu.tab.genus.n <- colnames(otu.tab.genus)[grepl("^NASH*",as.character(groups))] # After samples
otu.tab.genus.aldex <- data.frame(otu.tab.genus[,otu.tab.genus.h], otu.tab.genus[,otu.tab.genus.n]) # make a data frame
# make the vector of set membership in the same order as
otu.tab.genus.conds.aldex <- c(rep("Healthy", length(otu.tab.genus.h)), rep("NASH", length(otu.tab.genus.n)))
# generate 128 Dirichlet Monte-Carlo replicates
x.otu.tab.genus <- aldex.clr(otu.tab.genus.aldex, mc.samples=128, verbose=FALSE)
## [1] "operating in serial mode"
# calculate p values for each replicate and report the mean
x.otu.tab.genus.t <- aldex.ttest(x.otu.tab.genus, otu.tab.genus.conds.aldex)
# calculate mean effect sizes
x.otu.tab.genus.e <- aldex.effect(x.otu.tab.genus, otu.tab.genus.conds.aldex, verbose=FALSE)
## [1] "operating in serial mode"
# save it all in a data frame
x.otu.tab.genus.all <- data.frame(x.otu.tab.genus.e,x.otu.tab.genus.t)


pdf("aldex_plots.pdf")

layout(matrix(c(1,2,3,1,2,3),2,3, byrow=T), widths=c(5,2,2), height=c(4,4))
par(mar=c(5,4,4,1)+0.1)
aldex.plot(x.all, test="wilcox", cutoff=0.05, all.cex=0.8, called.cex=1)
plot(x.all$effect, x.all$wi.eBH, log="y", pch=19, main="Effect",
cex=0.5, xlab="Effect size", ylab="Expected Benjamini-Hochberg P")
abline(h=0.05, lty=2)
plot(x.all$diff.btw, x.all$wi.eBH, log="y", pch=19, main="Volcano",
cex=0.5, xlab="Difference", ylab="Expected Benjamini-Hochberg P")
abline(h=0.05, lty=2)

aldex.plot(x.filter.all, test="wilcox", cutoff=0.05, all.cex=0.8, called.cex=1)
plot(x.filter.all$effect, x.filter.all$wi.eBH, log="y", pch=19, main="Effect",
cex=0.5, xlab="Effect size", ylab="Expected Benjamini-Hochberg P")
abline(h=0.05, lty=2)
plot(x.filter.all$diff.btw, x.filter.all$wi.eBH, log="y", pch=19, main="Volcano",
cex=0.5, xlab="Difference", ylab="Expected Benjamini-Hochberg P")
abline(h=0.05, lty=2)

aldex.plot(x.genus.all, test="wilcox", cutoff=0.05, all.cex=0.8, called.cex=1)
plot(x.genus.all$effect, x.genus.all$wi.eBH, log="y", pch=19, main="Effect",
cex=0.5, xlab="Effect size", ylab="Expected Benjamini-Hochberg P")
abline(h=0.05, lty=2)
plot(x.genus.all$diff.btw, x.genus.all$wi.eBH, log="y", pch=19, main="Volcano",
cex=0.5, xlab="Difference", ylab="Expected Benjamini-Hochberg P")
abline(h=0.05, lty=2)

aldex.plot(x.otu.tab.genus.all, test="wilcox", cutoff=0.05, all.cex=0.8, called.cex=1)
plot(x.otu.tab.genus.all$effect, x.otu.tab.genus.all$wi.eBH, log="y", pch=19, main="Effect",
cex=0.5, xlab="Effect size", ylab="Expected Benjamini-Hochberg P")
abline(h=0.05, lty=2)
plot(x.otu.tab.genus.all$diff.btw, x.otu.tab.genus.all$wi.eBH, log="y", pch=19, main="Volcano",
cex=0.5, xlab="Difference", ylab="Expected Benjamini-Hochberg P")
abline(h=0.05, lty=2)

dev.off()

# COMPARE EFFECT SIZE WITH 16S


## sanity check to make sure all your counts have metadata
# which(!(colnames(otu.tab) %in% rownames(metadata)))

h.metnash <- otu.tab.genus
h.metnash.cond <- groups
h.metnash <- h.metnash[,which(h.metnash.cond == "Healthy Metagenomic" | h.metnash.cond == "NASH Metagenomic")]
h.metnash.cond <- h.metnash.cond[which(h.metnash.cond == "Healthy Metagenomic" | h.metnash.cond == "NASH Metagenomic")]

h.metnash.aldex <- aldex(data.frame(h.metnash),as.character(h.metnash.cond))


d.groups <- metadata$SSvsNASH[match(colnames(d.genus),rownames(metadata))]
d.originalgroups <- d.groups

d.cond <- d.groups
d.cond[which(is.na(d.cond))] <- "Healthy Metagenomic"
d.cond[which(d.cond==1)] <- "NASH Metagenomic"
d.cond <- as.factor(d.cond)

h.metnash.d <- d.genus

h.metnash.d.aldex <- aldex(data.frame(h.metnash.d),as.character(d.cond))

d.select <- d.genus
suspect.samples <- c("HLD_80","HLD_85","CL_165")
d.select <- d.genus[,which(!(colnames(d) %in% suspect.samples))]
d.select <- d.select[which(apply(d.select,1,sum)>0),]
d.select.cond <- colnames(d.select)
d.select.cond <- sub("CL.*$","NASH Metagenomic",d.select.cond)
d.select.cond <- sub("HLD.*$","Healthy Metagenomic",d.select.cond)

h.metnash.d.select.aldex <- aldex(data.frame(d.select),d.select.cond)

mycolor <- c(col2rgb("turquoise4"))
red <- mycolor[1]
green <- mycolor[2]
blue <- mycolor[3]
mycolor <- rgb(red/255, green/255, blue/255, 0.3)

pdf("metaphlan_vs_16S_effect_sizes.pdf")
common.genus <- rownames(h.metnash.d.aldex)[which(rownames(h.metnash.d.aldex) %in% rownames(h.metnash.aldex))]
otu.tab.common.effect <- h.metnash.aldex$effect[match(common.genus,rownames(h.metnash.aldex))]
d.common.effect <- h.metnash.d.aldex$effect[match(common.genus,rownames(h.metnash.d.aldex))]
plot(otu.tab.common.effect, d.common.effect, pch=19,col=mycolor, main="Effect sizes of healthy vs extreme NASH\nfor MetaPhlAn results vs. 16S sequencing",xlab="16S rRNA gene tag sequencing",ylab="MetaPhlAn results from metagenomic sequencing")
cor(otu.tab.common.effect, y = d.common.effect, use = "everything", method = "spearman")
# [1] 0.2136375

common.genus <- rownames(h.metnash.d.select.aldex)[which(rownames(h.metnash.d.select.aldex) %in% rownames(h.metnash.aldex))]
otu.tab.common.effect <- h.metnash.aldex$effect[match(common.genus,rownames(h.metnash.aldex))]
d.select.common.effect <- h.metnash.d.select.aldex$effect[match(common.genus,rownames(h.metnash.d.select.aldex))]
plot(otu.tab.common.effect, d.select.common.effect, pch=19,col=mycolor, main="Effect sizes of healthy vs extreme NASH\nfor MetaPhlAn results vs. 16S sequencing",xlab="16S rRNA gene tag sequencing",ylab="MetaPhlAn results from metagenomic sequencing")
cor(otu.tab.common.effect, y = d.select.common.effect, use = "everything", method = "spearman")
# [1] 0.268084
dev.off()

# see if 16S effect sizes correspond with qPCR
x.otu.tab.genus.all <- x.otu.tab.genus.all[order(abs(x.otu.tab.genus.all$effect),decreasing=TRUE),]
print(head(x.otu.tab.genus.all))

# rab.all rab.win.Healthy rab.win.NASH   diff.btw diff.win
# Lactobacillus  2.196815        1.333688     2.776964  1.9528688 2.657096
# Alistipes      3.485262        4.295601     2.592019 -1.1563853 2.054190
# Bacteroides    7.569359        7.934155     7.053593 -0.8011732 1.462026
# Blautia        8.483616        8.666417     8.238189 -0.4205494 0.786462
# Intestinimonas 3.025118        3.317382     2.776887 -0.6966185 1.410160
# Coprococcus    5.923939        6.289015     5.776760 -0.7112986 1.562239
# diff.btw.025 diff.btw.975 diff.win.025 diff.win.975     effect
# Lactobacillus     -6.105146    6.1835785    0.3979695     8.869206  0.6294914
# Alistipes         -4.816366    2.4503889    0.2638811     4.647072 -0.5801879
# Bacteroides       -2.970889    1.9698202    0.2478455     3.116608 -0.5111960
# Blautia           -1.860936    0.9820176    0.1605957     1.846667 -0.5068033
# Intestinimonas    -3.723969    1.7766566    0.2078378     3.508298 -0.4562087
# Coprococcus       -3.135996    1.9401785    0.2594317     3.173695 -0.4366168
# effect.025 effect.975   overlap       we.ep    we.eBH
# Lactobacillus   -2.595383   6.511377 0.2317708 0.007496508 0.1323528
# Alistipes       -6.258857   2.030277 0.2623698 0.002900021 0.1213582
# Bacteroides     -5.572032   2.007727 0.2864583 0.010575633 0.1394713
# Blautia         -4.950679   1.770154 0.2830189 0.007687336 0.1301380
# Intestinimonas  -5.840924   2.029911 0.2979831 0.006411170 0.1273803
# Coprococcus     -4.611679   2.070917 0.3350683 0.017021931 0.1683449
#       wi.ep     wi.eBH
# Lactobacillus  0.0009600779 0.06219896
# Alistipes      0.0031407310 0.09884836
# Bacteroides    0.0103019133 0.16021214
# Blautia        0.0114378306 0.15939895
# Intestinimonas 0.0149239237 0.18027105
# Coprococcus    0.0420344247 0.27472029

h.metnash.aldex <- h.metnash.aldex[order(abs(h.metnash.aldex$effect),decreasing=TRUE),]
print(head(h.metnash.aldex))

# rab.all rab.win.Healthy.Metagenomic
# Ruminococcus           6.8726076                   7.5502488
# Paraprevotella        -0.1205997                  -1.4204805
# Phascolarctobacterium  3.0790176                   1.3912225
# Alistipes              3.2401873                   3.8960106
# Coprococcus            6.1074078                   6.4200064
# Odoribacter           -1.1875189                   0.1572045
# rab.win.NASH.Metagenomic   diff.btw diff.win diff.btw.025
# Ruminococcus                         6.2909093 -1.1668231 1.188396    -2.981244
# Paraprevotella                       0.9295206  2.6854951 3.216184    -2.104623
# Phascolarctobacterium                4.4343099  1.8423307 2.381068    -2.409790
# Alistipes                            2.4842489 -1.2048407 1.888123    -4.368788
# Coprococcus                          5.9638318 -0.8262721 1.253413    -2.915490
# Odoribacter                         -2.0470168 -1.9224693 3.220999   -10.697790
# diff.btw.975 diff.win.025 diff.win.975     effect
# Ruminococcus             1.3771325    0.1851817     2.708947 -0.8798468
# Paraprevotella          10.8807941    0.3907753    10.432557  0.8664093
# Phascolarctobacterium    5.2531514    0.2888974     4.292936  0.8585089
# Alistipes                1.7327370    0.2702767     3.652876 -0.6780909
# Coprococcus              0.9949517    0.1482030     2.254421 -0.6318341
# Odoribacter              3.5691321    0.5819410    10.986667 -0.5911948
# effect.025 effect.975   overlap      we.ep    we.eBH
# Ruminococcus          -7.4892279   1.948368 0.2043683 0.01381314 0.4060738
# Paraprevotella        -0.8612033  10.134242 0.1359377 0.01112611 0.3762764
# Phascolarctobacterium -1.7039929   9.530847 0.1843751 0.01538188 0.3996781
# Alistipes             -6.2947490   1.401327 0.2046876 0.02567951 0.4512520
# Coprococcus           -7.9097586   1.663411 0.2683308 0.03642488 0.5078688
# Odoribacter           -7.6168979   1.443395 0.2218751 0.04875179 0.5538504
#    wi.ep    wi.eBH
# Ruminococcus          0.022090876 0.4283251
# Paraprevotella        0.005509462 0.3103027
# Phascolarctobacterium 0.019750965 0.4262033
# Alistipes             0.023957104 0.4411569
# Coprococcus           0.091530772 0.7174758
# Odoribacter           0.043399808 0.5202808


### TOP AND BOTTOM DECILE OF OTUs
top <- c(36, 84, 120, 538, 635, 422, 15, 88, 174, 318, 178, 142, 1000, 937, 64, 65, 211, 1447, 1279, 742, 280, 1245, 332, 1331, 286, 155, 170, 185, 878, 116, 100, 197, 30, 338, 160, 66, 297, 350, 225, 293, 96, 387, 631, 112, 1097, 261, 72, 236, 359, 61, 758, 990, 80)
bottom <- c(28, 10, 316, 57, 254, 344, 7, 127, 5, 1289, 1118, 259, 257, 83, 103, 50, 91, 846, 329, 202, 122, 287, 70, 195, 110, 558, 157, 17, 383, 1101, 815, 1005, 79, 996, 180, 40, 0, 23, 54, 16, 2, 145, 224, 14, 320, 8, 1202, 18, 111, 12, 718, 31, 176)
top <- as.character(top)
bottom <- as.character(bottom)

top.genus <- otu.genus[match(top,rownames(otu.tab))]
print(top.genus)
# [1] "Phascolarctobacterium" "Lactobacillus"         "Paraprevotella"       
# [4] "Incertae_Sedis"        "Marvinbryantia"        "Incertae_Sedis"       
# [7] "Bifidobacterium"       "Incertae_Sedis"        "Paraprevotella"       
# [10] "unclassified"          "unclassified"          "Butyricicoccus"       
# [13] "Incertae_Sedis"        "Ruminococcus"          "unclassified"         
# [16] "unclassified"          "Lactobacillus"         "Subdoligranulum"      
# [19] "Incertae_Sedis"        "unclassified"          "Olsenella"            
# [22] "Subdoligranulum"       "unclassified"          "Prevotella"           
# [25] "Anaerostipes"          "unclassified"          "unclassified"         
# [28] "Incertae_Sedis"        "Faecalibacterium"      "Sutterella"           
# [31] "unclassified"          "Alistipes"             "Roseburia"            
# [34] "unclassified"          "Enterorhabdus"         "Ruminococcus"         
# [37] "unclassified"          "Dialister"             "unclassified"         
# [40] "Incertae_Sedis"        "Bacteroides"           "unclassified"         
# [43] "unclassified"          "Desulfovibrio"         "unclassified"         
# [46] "Blautia"               "Blautia"               "unclassified"         
# [49] "unclassified"          "Parasutterella"        "Incertae_Sedis"       
# [52] "Subdoligranulum"       "Acidaminococcus"     

bottom.genus <- otu.genus[match(bottom,rownames(otu.tab))]
print(bottom.genus)
# [1] "Incertae_Sedis"     "Akkermansia"        "Parabacteroides"   
# [4] "Alistipes"          "Incertae_Sedis"     "unclassified"      
# [7] "Streptococcus"      "unclassified"       "Dorea"             
# [10] "Roseburia"          "Incertae_Sedis"     "unclassified"      
# [13] "unclassified"       "Blautia"            "unclassified"      
# [16] "Roseburia"          "Odoribacter"        "Incertae_Sedis"    
# [19] "unclassified"       "Odoribacter"        "Turicibacter"      
# [22] "unclassified"       "Ruminococcus"       "unclassified"      
# [25] "unclassified"       "Bacteroides"        "Incertae_Sedis"    
# [28] "Dialister"          "Bacteroides"        "Alloprevotella"    
# [31] "Bacteroides"        "Alistipes"          "Adlercreutzia"     
# [34] "unclassified"       "Anaerovorax"        "unclassified"      
# [37] "Pseudobutyrivibrio" "Bacteroides"        "Bacteroides"       
# [40] "Blautia"            "Faecalibacterium"   "unclassified"      
# [43] "Incertae_Sedis"     "Subdoligranulum"    "Incertae_Sedis"    
# [46] "Ruminococcus"       "Incertae_Sedis"     "Coprococcus"       
# [49] "unclassified"       "unclassified"       "Subdoligranulum"   
# [52] "Dorea"              "Incertae_Sedis"       

both.genus <- top.genus[which(top.genus %in% bottom.genus)]
top.genus.only <- top.genus[which(! top.genus %in% bottom.genus)]
bottom.genus.only <- bottom.genus[which(! bottom.genus %in% top.genus)]

pdf("metaphlan_vs_16S_effect_sizes_decile_genus.pdf")
common.genus <- rownames(h.metnash.d.aldex)[which(rownames(h.metnash.d.aldex) %in% rownames(h.metnash.aldex))]
genus.groups <- rep("normal",length(common.genus))

print(top.genus[which(! top.genus %in% common.genus)])
# [1] "Incertae_Sedis" "Incertae_Sedis" "Incertae_Sedis" "unclassified"  
# [5] "unclassified"   "Incertae_Sedis" "unclassified"   "unclassified"  
# [9] "Incertae_Sedis" "unclassified"   "unclassified"   "unclassified"  
# [13] "unclassified"   "Incertae_Sedis" "unclassified"   "unclassified"  
# [17] "unclassified"   "unclassified"   "Incertae_Sedis" "unclassified"  
# [21] "unclassified"   "unclassified"   "unclassified"   "unclassified"  
# [25] "Incertae_Sedis"
print(bottom.genus[which(! bottom.genus %in% common.genus)])
# [1] "Incertae_Sedis"     "Incertae_Sedis"     "unclassified"      
# [4] "unclassified"       "Incertae_Sedis"     "unclassified"      
# [7] "unclassified"       "unclassified"       "Incertae_Sedis"    
# [10] "unclassified"       "unclassified"       "unclassified"      
# [13] "unclassified"       "Incertae_Sedis"     "Alloprevotella"    
# [16] "unclassified"       "Anaerovorax"        "unclassified"      
# [19] "Pseudobutyrivibrio" "unclassified"       "Incertae_Sedis"    
# [22] "Incertae_Sedis"     "Incertae_Sedis"     "unclassified"      
# [25] "unclassified"       "Incertae_Sedis"   
genus.groups[which(common.genus %in% both.genus)] <- "both"
genus.groups[which(common.genus %in% top.genus.only)] <- "top"
genus.groups[which(common.genus %in% bottom.genus.only)] <- "bottom"
genus.groups <- as.factor(genus.groups)

# make plots with top and bottom decile effect sizes from metnash comparison colored
firebrick <- c(col2rgb("firebrick"))
red <- firebrick[1]
green <- firebrick[2]
blue <- firebrick[3]
firebrick <- rgb(red/255, green/255, blue/255, 0.3)

violet <- c(col2rgb("darkviolet"))
red <- violet[1]
green <- violet[2]
blue <- violet[3]
violet <- rgb(red/255, green/255, blue/255, 0.3)

myblue <- c(col2rgb("deepskyblue3"))
red <- myblue[1]
green <- myblue[2]
blue <- myblue[3]
myblue <- rgb(red/255, green/255, blue/255, 0.3)

mygray <- c(col2rgb("gray18"))
red <- mygray[1]
green <- mygray[2]
blue <- mygray[3]
mygray <- rgb(red/255, green/255, blue/255, 0.3)

palette(c(violet,myblue,mygray,firebrick))

otu.tab.common.effect <- h.metnash.aldex$effect[match(common.genus,rownames(h.metnash.aldex))]
d.common.effect <- h.metnash.d.aldex$effect[match(common.genus,rownames(h.metnash.d.aldex))]
plot(otu.tab.common.effect, d.common.effect, pch=19,col=genus.groups, main="Effect sizes of healthy vs extreme NASH\nfor MetaPhlAn results vs. 16S sequencing",xlab="16S rRNA gene tag sequencing",ylab="MetaPhlAn results from metagenomic sequencing")
cor(otu.tab.common.effect, y = d.common.effect, use = "everything", method = "spearman")
# [1] 0.2021001

common.genus <- rownames(h.metnash.d.select.aldex)[which(rownames(h.metnash.d.select.aldex) %in% rownames(h.metnash.aldex))]
otu.tab.common.effect <- h.metnash.aldex$effect[match(common.genus,rownames(h.metnash.aldex))]
d.select.common.effect <- h.metnash.d.select.aldex$effect[match(common.genus,rownames(h.metnash.d.select.aldex))]
plot(otu.tab.common.effect, d.select.common.effect, pch=19,col=genus.groups, main="Effect sizes of healthy vs extreme NASH\nfor MetaPhlAn results vs. 16S sequencing",xlab="16S rRNA gene tag sequencing",ylab="MetaPhlAn results from metagenomic sequencing")
cor(otu.tab.common.effect, y = d.select.common.effect, use = "everything", method = "spearman")
# [1] 0.235092
dev.off()

print(paste(round(h.metnash.aldex$effect[match(top.genus, rownames(h.metnash.aldex))],digits=3),collapse="\n"))
# [1] "0.827\n0.587\n0.843\nNA\n0.077\nNA\n0.188\nNA\n0.843\nNA\nNA\n0.449\nNA\n-0.866\nNA\nNA\n0.587\n-0.087\nNA\nNA\n0.318\n-0.087\nNA\n0.212\n0.244\nNA\nNA\nNA\n-0.226\n0.162\nNA\n-0.687\n0.18\nNA\n0.429\n-0.866\nNA\n-0.297\nNA\nNA\n-0.356\nNA\nNA\n0.283\nNA\n-0.031\n-0.031\nNA\nNA\n0.389\nNA\n-0.087\n0.345"

print(paste(round(h.metnash.aldex$effect[match(bottom.genus, rownames(h.metnash.aldex))],digits=3),collapse="\n"))
# [1] "NA\n-0.496\n-0.313\n-0.687\nNA\nNA\n-0.18\nNA\n-0.267\n0.18\nNA\nNA\nNA\n-0.031\nNA\n0.18\n-0.541\nNA\nNA\n-0.541\n-0.52\nNA\n-0.866\nNA\nNA\n-0.356\nNA\n-0.297\n-0.356\n-0.001\n-0.356\n-0.687\n-0.372\nNA\n-0.476\nNA\n-0.423\n-0.356\n-0.356\n-0.031\n-0.226\nNA\nNA\n-0.087\nNA\n-0.866\nNA\n-0.647\nNA\nNA\n-0.087\n-0.267\nNA"

print(paste(round(h.metnash.d.aldex$effect[match(top.genus, rownames(h.metnash.aldex))],digits=3),collapse="\n"))
# [1] "0.081\n-0.899\n0.064\nNA\n0.269\nNA\n0.032\nNA\n0.064\nNA\nNA\n0.28\nNA\n0.023\nNA\nNA\n-0.899\n-0.177\nNA\nNA\n0.141\n-0.177\nNA\n0.188\n0.344\nNA\nNA\nNA\n-0.173\n0.051\nNA\n0.053\n0.168\nNA\n0.224\n0.023\nNA\n-0.038\nNA\nNA\n-0.124\nNA\nNA\n0.046\nNA\n0.192\n0.192\nNA\nNA\n0.521\nNA\n-0.177\n-0.737"

print(paste(round(h.metnash.d.aldex$effect[match(bottom.genus, rownames(h.metnash.aldex))],digits=3),collapse="\n"))
# [1] "NA\n0.343\n0.016\n0.053\nNA\nNA\n0.233\nNA\n-0.154\n0.168\nNA\nNA\nNA\n0.192\nNA\n0.168\n-0.333\nNA\nNA\n-0.333\n-0.717\nNA\n0.023\nNA\nNA\n-0.124\nNA\n-0.038\n-0.124\n0.863\n-0.124\n0.053\n-0.298\nNA\n-0.152\nNA\n-0.091\n-0.124\n-0.124\n0.192\n-0.173\nNA\nNA\n-0.177\nNA\n0.023\nNA\n0.469\nNA\nNA\n-0.177\n-0.154\nNA"

top.taxonomy <- taxonomy[match(top,rownames(otu.tab))]
print(top.taxonomy)
# [1] Bacteria;Firmicutes;Negativicutes;Selenomonadales;Acidaminococcaceae;Phascolarctobacterium;|100      
# [2] Bacteria;Firmicutes;Bacilli;Lactobacillales;Lactobacillaceae;Lactobacillus;|97                       
# [3] Bacteria;Bacteroidetes;Bacteroidia;Bacteroidales;Prevotellaceae;Paraprevotella;|100                  
# [4] Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae;Incertae_Sedis;|98                      
# [5] Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae;Marvinbryantia;|77                      
# [6] Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae;Incertae_Sedis;|73                      
# [7] Bacteria;Actinobacteria;Actinobacteria;Bifidobacteriales;Bifidobacteriaceae;Bifidobacterium;|100     
# [8] Bacteria;Firmicutes;Clostridia;Clostridiales;Ruminococcaceae;Incertae_Sedis;|72                      
# [9] Bacteria;Bacteroidetes;Bacteroidia;Bacteroidales;Prevotellaceae;Paraprevotella;|100                  
# [10] Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae;unclassified;|100                       
# [11] Bacteria;Firmicutes;Clostridia;Clostridiales;unclassified;unclassified;|72                           
# [12] Bacteria;Firmicutes;Clostridia;Clostridiales;Ruminococcaceae;Butyricicoccus;|71                      
# [13] Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae;Incertae_Sedis;|91                      
# [14] Bacteria;Firmicutes;Clostridia;Clostridiales;Ruminococcaceae;Ruminococcus;|93                        
# [15] Bacteria;Actinobacteria;Coriobacteriia;Coriobacteriales;Coriobacteriaceae;unclassified;|97           
# [16] Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae;unclassified;|98                        
# [17] Bacteria;Firmicutes;Bacilli;Lactobacillales;Lactobacillaceae;Lactobacillus;|98                       
# [18] Bacteria;Firmicutes;Clostridia;Clostridiales;Ruminococcaceae;Subdoligranulum;|87                     
# [19] Bacteria;Firmicutes;Clostridia;Clostridiales;Ruminococcaceae;Incertae_Sedis;|98                      
# [20] Bacteria;Firmicutes;Clostridia;Clostridiales;Ruminococcaceae;unclassified;|100                       
# [21] Bacteria;Actinobacteria;Coriobacteriia;Coriobacteriales;Coriobacteriaceae;Olsenella;|91              
# [22] Bacteria;Firmicutes;Clostridia;Clostridiales;Ruminococcaceae;Subdoligranulum;|98                     
# [23] Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae;unclassified;|100                       
# [24] Bacteria;Bacteroidetes;Bacteroidia;Bacteroidales;Prevotellaceae;Prevotella;|99                       
# [25] Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae;Anaerostipes;|100                       
# [26] Bacteria;Bacteroidetes;Bacteroidia;Bacteroidales;Prevotellaceae;unclassified;|70                     
# [27] Bacteria;Firmicutes;Clostridia;Clostridiales;unclassified;unclassified;|92                           
# [28] Bacteria;Firmicutes;Clostridia;Clostridiales;Ruminococcaceae;Incertae_Sedis;|99                      
# [29] Bacteria;Firmicutes;Clostridia;Clostridiales;Ruminococcaceae;Faecalibacterium;|100                   
# [30] Bacteria;Proteobacteria;Betaproteobacteria;Burkholderiales;Alcaligenaceae;Sutterella;|100            
# [31] Bacteria;Firmicutes;Clostridia;Clostridiales;unclassified;unclassified;|73                           
# [32] Bacteria;Bacteroidetes;Bacteroidia;Bacteroidales;Rikenellaceae;Alistipes;|100                        
# [33] Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae;Roseburia;|98                           
# [34] Bacteria;Bacteroidetes;Bacteroidia;Bacteroidales;Prevotellaceae;unclassified;|75                     
# [35] Bacteria;Actinobacteria;Coriobacteriia;Coriobacteriales;Coriobacteriaceae;Enterorhabdus;|72          
# [36] Bacteria;Firmicutes;Clostridia;Clostridiales;Ruminococcaceae;Ruminococcus;|100                       
# [37] Bacteria;Firmicutes;Clostridia;Clostridiales;unclassified;unclassified;|98                           
# [38] Bacteria;Firmicutes;Negativicutes;Selenomonadales;Veillonellaceae;Dialister;|100                     
# [39] Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae;unclassified;|100                       
# [40] Bacteria;Firmicutes;Clostridia;Clostridiales;Family_XIII;Incertae_Sedis;|100                         
# [41] Bacteria;Bacteroidetes;Bacteroidia;Bacteroidales;Bacteroidaceae;Bacteroides;|100                     
# [42] Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae;unclassified;|100                       
# [43] Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae;unclassified;|100                       
# [44] Bacteria;Proteobacteria;Deltaproteobacteria;Desulfovibrionales;Desulfovibrionaceae;Desulfovibrio;|100
# [45] Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae;unclassified;|100                       
# [46] Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae;Blautia;|96                             
# [47] Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae;Blautia;|97                             
# [48] Bacteria;Firmicutes;Clostridia;Clostridiales;Ruminococcaceae;unclassified;|100                       
# [49] Bacteria;Firmicutes;Clostridia;Clostridiales;Christensenellaceae;unclassified;|99                    
# [50] Bacteria;Proteobacteria;Betaproteobacteria;Burkholderiales;Alcaligenaceae;Parasutterella;|100        
# [51] Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae;Incertae_Sedis;|100                     
# [52] Bacteria;Firmicutes;Clostridia;Clostridiales;Ruminococcaceae;Subdoligranulum;|92                     
# [53] Bacteria;Firmicutes;Negativicutes;Selenomonadales;Acidaminococcaceae;Acidaminococcus;|100   
    
bottom.taxonomy <- taxonomy[match(bottom,rownames(otu.tab))]
print(bottom.taxonomy)
# [1] Bacteria;Firmicutes;Clostridia;Clostridiales;Ruminococcaceae;Incertae_Sedis;|100                 
# [2] Bacteria;Verrucomicrobia;Verrucomicrobiae;Verrucomicrobiales;Verrucomicrobiaceae;Akkermansia;|100
# [3] Bacteria;Bacteroidetes;Bacteroidia;Bacteroidales;Porphyromonadaceae;Parabacteroides;|100         
# [4] Bacteria;Bacteroidetes;Bacteroidia;Bacteroidales;Rikenellaceae;Alistipes;|100                    
# [5] Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae;Incertae_Sedis;|73                  
# [6] Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae;unclassified;|100                   
# [7] Bacteria;Firmicutes;Bacilli;Lactobacillales;Streptococcaceae;Streptococcus;|100                  
# [8] Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae;unclassified;|100                   
# [9] Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae;Dorea;|100                          
# [10] Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae;Roseburia;|83                       
# [11] Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae;Incertae_Sedis;|82                  
# [12] Bacteria;Firmicutes;Clostridia;Clostridiales;Christensenellaceae;unclassified;|98                
# [13] Bacteria;Firmicutes;Clostridia;Clostridiales;Christensenellaceae;unclassified;|100               
# [14] Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae;Blautia;|93                         
# [15] Bacteria;Firmicutes;Clostridia;Clostridiales;Ruminococcaceae;unclassified;|100                   
# [16] Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae;Roseburia;|91                       
# [17] Bacteria;Bacteroidetes;Bacteroidia;Bacteroidales;Porphyromonadaceae;Odoribacter;|100             
# [18] Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae;Incertae_Sedis;|85                  
# [19] Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae;unclassified;|98                    
# [20] Bacteria;Bacteroidetes;Bacteroidia;Bacteroidales;Porphyromonadaceae;Odoribacter;|100             
# [21] Bacteria;Firmicutes;Erysipelotrichia;Erysipelotrichales;Erysipelotrichaceae;Turicibacter;|100    
# [22] Bacteria;Firmicutes;Clostridia;Clostridiales;Christensenellaceae;unclassified;|100               
# [23] Bacteria;Firmicutes;Clostridia;Clostridiales;Ruminococcaceae;Ruminococcus;|100                   
# [24] Bacteria;Firmicutes;Clostridia;Clostridiales;Christensenellaceae;unclassified;|99                
# [25] Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae;unclassified;|100                   
# [26] Bacteria;Bacteroidetes;Bacteroidia;Bacteroidales;Bacteroidaceae;Bacteroides;|100                 
# [27] Bacteria;Firmicutes;Clostridia;Clostridiales;Ruminococcaceae;Incertae_Sedis;|84                  
# [28] Bacteria;Firmicutes;Negativicutes;Selenomonadales;Veillonellaceae;Dialister;|100                 
# [29] Bacteria;Bacteroidetes;Bacteroidia;Bacteroidales;Bacteroidaceae;Bacteroides;|100                 
# [30] Bacteria;Bacteroidetes;Bacteroidia;Bacteroidales;Prevotellaceae;Alloprevotella;|100              
# [31] Bacteria;Bacteroidetes;Bacteroidia;Bacteroidales;Bacteroidaceae;Bacteroides;|100                 
# [32] Bacteria;Bacteroidetes;Bacteroidia;Bacteroidales;Rikenellaceae;Alistipes;|100                    
# [33] Bacteria;Actinobacteria;Coriobacteriia;Coriobacteriales;Coriobacteriaceae;Adlercreutzia;|100     
# [34] Bacteria;Firmicutes;Clostridia;Clostridiales;Ruminococcaceae;unclassified;|100                   
# [35] Bacteria;Firmicutes;Clostridia;Clostridiales;Family_XIII;Anaerovorax;|91                         
# [36] Bacteria;Firmicutes;Clostridia;Clostridiales;Ruminococcaceae;unclassified;|76                    
# [37] Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae;Pseudobutyrivibrio;|98              
# [38] Bacteria;Bacteroidetes;Bacteroidia;Bacteroidales;Bacteroidaceae;Bacteroides;|100                 
# [39] Bacteria;Bacteroidetes;Bacteroidia;Bacteroidales;Bacteroidaceae;Bacteroides;|100                 
# [40] Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae;Blautia;|85                         
# [41] Bacteria;Firmicutes;Clostridia;Clostridiales;Ruminococcaceae;Faecalibacterium;|100               
# [42] Bacteria;Firmicutes;Clostridia;Clostridiales;unclassified;unclassified;|85                       
# [43] Bacteria;Firmicutes;Clostridia;Clostridiales;Ruminococcaceae;Incertae_Sedis;|100                 
# [44] Bacteria;Firmicutes;Clostridia;Clostridiales;Ruminococcaceae;Subdoligranulum;|99                 
# [45] Bacteria;Firmicutes;Clostridia;Clostridiales;Ruminococcaceae;Incertae_Sedis;|97                  
# [46] Bacteria;Firmicutes;Clostridia;Clostridiales;Ruminococcaceae;Ruminococcus;|100                   
# [47] Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae;Incertae_Sedis;|85                  
# [48] Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae;Coprococcus;|91                     
# [49] Bacteria;Firmicutes;Clostridia;Clostridiales;Ruminococcaceae;unclassified;|74                    
# [50] Bacteria;Firmicutes;Erysipelotrichia;Erysipelotrichales;Erysipelotrichaceae;unclassified;|100    
# [51] Bacteria;Firmicutes;Clostridia;Clostridiales;Ruminococcaceae;Subdoligranulum;|85                 
# [52] Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae;Dorea;|100                          
# [53] Bacteria;Firmicutes;Clostridia;Clostridiales;Family_XIII;Incertae_Sedis;|91  