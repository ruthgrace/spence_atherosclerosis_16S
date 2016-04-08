library(zCompositions)
library(randomcoloR)
library(compositions)
library(ALDEx2)
library(stringr)

# read metadata for 16S samples
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

colnames(metadata) <- gsub(" ",".",colnames(metadata))

# get 16S rRNA gene sequencing count
otu.tab <- read.table("data/td_OTU_tag_mapped_lineage.txt", header=T, row.names=1, sep="\t",quote="",comment.char="",stringsAsFactors=FALSE,check.names=FALSE)

# make samples consistent between metadata and count data
metadata.samples <- match(colnames(otu.tab),rownames(metadata))
metadata.samples <- metadata.samples[which(!is.na(metadata.samples))]
metadata <- metadata[metadata.samples,]

taxonomy <- otu.tab$taxonomy

otu.tab <- otu.tab[,c(1:(ncol(otu.tab)-1))]
otu.tab <- otu.tab[,which(colnames(otu.tab) %in% rownames(metadata))]

otu.tab <- as.matrix(t(otu.tab))

# remove all features with less than 1% abundance in all samples
otu.tab.sum <- apply(otu.tab,1,sum)
otu.tab.sum.threshhold <- otu.tab.sum*0.01
otu.tab.filter <- apply(otu.tab,2, function(x) { return(length(which(x > otu.tab.sum.threshhold))) } )
otu.tab <- otu.tab[,which(otu.tab.filter > 0)]

taxonomy <- taxonomy[which(otu.tab.filter > 0)]

otu.sum <- apply(otu.tab,2,sum)
otu.tab <- otu.tab[,order(otu.sum,decreasing=TRUE)]
taxonomy <- taxonomy[order(otu.sum,decreasing=TRUE)]

original.data <- otu.tab

otu.sp <- c(as.character(taxonomy))

for (i in c(1:length(taxonomy))) {
  otu.sp[i] <- paste(strsplit(otu.sp[i],c(";"))[[1]][c(6,7)],collapse="")
}

residual.order <- order(metadata$Standardized.Residual,decreasing=TRUE)
decile <- round(nrow(metadata)/10)
top.decile <- residual.order[c(1:decile)]
bottom.decile <- residual.order[c(nrow(metadata)-decile+1):nrow(metadata)]
otu.tab.extreme <- otu.tab[c(top.decile,bottom.decile),]

# conditions: Originally 0 meant steatohepatosis, and 1 meant NASH
residuals <- metadata$Standardized.Residual
groups <- rep("NA",length(residuals))
groups[which(residuals < -1)] <- "Protected"
groups[which(residuals >= -1 & residuals < 1)] <- "Explained"
groups[which(residuals >= 1)] <- "Unexplained"
groups <- as.factor(groups)

groups.extreme <- c(rep("Top decile",decile),rep("Bottom decile",decile))
groups.extreme <- as.factor(groups.extreme)

d <- t(otu.tab)
d.extreme <- t(otu.tab.extreme)
rownames(d) <- otu.sp
rownames(d.extreme) <- otu.sp
# it just happens that there are no zero count features in d.extreme, so it's not necessary to filter

# adjust zeros
d.adj.zero <- t(cmultRepl(t(d),method="CZM"))
d.extreme.adj.zero <- t(cmultRepl(t(d.extreme),method="CZM"))

d.names <- rownames(d.adj.zero)
d.extreme.names <- rownames(d.extreme.adj.zero)

taxa.col <- data.frame(as.character(rownames(d)),rownames(d))
colnames(taxa.col) <- c("taxon","color")
taxa.col[,2] <- distinctColorPalette(length(taxa.col[,2]))

conds <- data.frame(as.character(groups))
colnames(conds) <- "cond"

palette=palette(c(rgb(1,0,0,0.6), "green","black", rgb(0,1,1,0.6)))

my.col <- rep("NA",length(groups))
my.col[which(groups == "Explained")] <- "purple"
my.col[which(groups == "Protected")] <- "blue"
my.col[which(groups == "Unexplained")] <- "red"

my.extreme.col <- c(rep("red",20),rep("blue",20))

plot.biplot <- function(d.adj.zero,my.col,groups,legend.col) {
  d.prop <- apply(d.adj.zero,2,function(x){x/sum(x)})
  d.clr <- t(apply(d.prop,2,function(x){log(x) - mean(log(x))}))
  d.pcx <- prcomp(d.clr)
  layout(matrix(c(1,2),1,2, byrow=T), widths=c(6,2), heights=c(8,3))
  par(mgp=c(2,0.5,0), xpd=TRUE)
  # make a covariance biplot of the data with compositions function
  coloredBiplot(d.pcx, cex=c(0.6, 0.6),
  col=rgb(0,0,0,0.5),
  arrow.len=0.05,
  xlab=paste("PC1 ", round (sum(d.pcx$sdev[1]^2)/mvar(d.clr),3), sep=""),
  ylab=paste("PC2 ", round (sum(d.pcx$sdev[2]^2)/mvar(d.clr),3), sep=""),
  xlabs.col=my.col,
  expand=0.8,var.axes=FALSE, scale=1, main="Biplot")
  legend("bottomleft",inset=c(-0.15,-0.25),levels(groups),pch=19,col=legend.col)
  barplot(d.pcx$sdev^2/mvar(d.clr),  ylab="variance explained", xlab="Component", main="Scree plot")
  return(d.clr)
}

pdf("biplots.pdf")

d.clr <- plot.biplot(d.adj.zero,my.col,groups,c("purple","blue","red"))

d.extreme.clr <- plot.biplot(d.extreme.adj.zero,my.extreme.col,groups.extreme,c("blue","red"))

dev.off()

plot.dendogram.barplot <- function(d.clr, taxa.col, d.names, groups) {
  d.dist <- dist(d.clr, method="euclidian")
  d.conditions <- as.character(groups)
  d.conditions <- gsub(" ","_",d.conditions)
  attributes(d.dist)$Labels <- paste(d.conditions, attributes(d.dist)$Labels, sep="_")
  d.hc <- hclust(d.dist, method="ward.D2")
  d.order <- d.adj.zero[,d.hc$order]
  d.acomp <- acomp(t(d.order))

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
}

pdf("dendogram_barplot.pdf")

plot.dendogram.barplot(d.clr,taxa.col,d.names,groups)

plot.dendogram.barplot(d.extreme.clr,taxa.col,d.extreme.names,groups.extreme)

dev.off()

plot.aldex.volcano <- function(d.aldex,d.conds) {
  x <- aldex.clr(d.aldex, mc.samples=128, verbose=FALSE)
  ## [1] "operating in serial mode"
  # calculate p values for each replicate and report the mean
  x.t <- aldex.ttest(x, d.conds)
  # calculate mean effect sizes
  x.e <- aldex.effect(x, d.conds, verbose=FALSE)
  ## [1] "operating in serial mode"
  # save it all in a data frame
  x.all <- data.frame(x.e,x.t)
  
  layout(matrix(c(1,2,3,1,2,3),2,3, byrow=T), widths=c(5,2,2), height=c(4,4))
  par(mar=c(5,4,4,1)+0.1)
  aldex.plot(x.all, test="wilcox", cutoff=0.05, all.cex=0.8, called.cex=1)
  plot(x.all$effect, x.all$wi.eBH, log="y", pch=19, main="Effect",
  cex=0.5, xlab="Effect size", ylab="Expected Benjamini-Hochberg P")
  abline(h=0.05, lty=2)
  plot(x.all$diff.btw, x.all$wi.eBH, log="y", pch=19, main="Volcano",
  cex=0.5, xlab="Difference", ylab="Expected Benjamini-Hochberg P")
  abline(h=0.05, lty=2)
}





# generate the dataset by making a data frame of
ep.h <- colnames(d)[which(groups == "Explained")] # Before samples
ep.n <- colnames(d)[which(groups == "Protected")] # After samples
ep.aldex <- data.frame(d[,ep.h], d[,ep.n],check.names=FALSE) # make a data frame
# make the vector of set membership in the same order as
ep.conds.aldex <- c(rep("Explained", length(ep.h)), rep("Protected", length(ep.n)))

# generate the dataset by making a data frame of
eu.h <- colnames(d)[which(groups == "Explained")] # Before samples
eu.n <- colnames(d)[which(groups == "Unexplained")] # After samples
eu.aldex <- data.frame(d[,eu.h], d[,eu.n],check.names=FALSE) # make a data frame
# make the vector of set membership in the same order as
eu.conds.aldex <- c(rep("Explained", length(eu.h)), rep("Unexplained", length(eu.n)))

# generate the dataset by making a data frame of
up.h <- colnames(d)[which(groups == "Unexplained")] # Before samples
up.n <- colnames(d)[which(groups == "Protected")] # After samples
up.aldex <- data.frame(d[,up.h], d[,up.n],check.names=FALSE) # make a data frame
# make the vector of set membership in the same order as
up.conds.aldex <- c(rep("Unexplained", length(up.h)), rep("Protected", length(up.n)))

# generate the dataset by making a data frame of
d.extreme.aldex <- data.frame(d.extreme) # make a data frame
# make the vector of set membership in the same order as
extreme.conds.aldex <- as.character(groups.extreme)

pdf("aldex_plots.pdf")

plot.aldex.volcano(ep.aldex,ep.conds.aldex)
plot.aldex.volcano(eu.aldex,eu.conds.aldex)
plot.aldex.volcano(up.aldex,up.conds.aldex)
plot.aldex.volcano(d.extreme.aldex,extreme.conds.aldex)

dev.off()

## sanity check to make sure all your counts have metadata
# which(!(colnames(otu.tab) %in% rownames(metadata)))

h.ep.aldex <- aldex(data.frame(ep.aldex),as.character(ep.conds.aldex))
h.eu.aldex <- aldex(data.frame(eu.aldex),as.character(eu.conds.aldex))
h.up.aldex <- aldex(data.frame(up.aldex),as.character(up.conds.aldex))
h.extreme.aldex <- aldex(data.frame(d.extreme.aldex),as.character(extreme.conds.aldex))

mycolor <- c(col2rgb("turquoise4"))
red <- mycolor[1]
green <- mycolor[2]
blue <- mycolor[3]
mycolor <- rgb(red/255, green/255, blue/255, 0.3)

pdf("16S_effect_sizes.pdf")
plot(h.ep.aldex$effect, h.up.aldex$effect, pch=19,col=mycolor, main="Effect sizes of Explained vs. Protected\nand Unexplained vs. Protected",xlab="Explained vs. Protected",ylab="Unexplained vs. Protected")
cor(h.ep.aldex$effect, y = h.up.aldex$effect, use = "everything", method = "spearman")

plot(h.eu.aldex$effect, h.up.aldex$effect, pch=19,col=mycolor, main="Effect sizes of Explained vs. Unexplained\nand Unexplained vs. Protected",xlab="Explained vs. Unexplained",ylab="Unexplained vs. Protected")
cor(h.eu.aldex$effect, y = h.up.aldex$effect, use = "everything", method = "spearman")

plot(h.up.aldex$effect, h.extreme.aldex$effect, pch=19,col=mycolor, main="Effect sizes of Unexplained vs. Protected\nand Top decile vs. Bottom decile subset",xlab="Unexplained vs. Protected",ylab="Top decile vs. Bottom decile subset")
cor(h.up.aldex$effect, y = h.extreme.aldex$effect, use = "everything", method = "spearman")
dev.off()

h.ep.aldex <- h.ep.aldex[order(abs(h.ep.aldex$effect),decreasing=TRUE),]
h.eu.aldex <- h.eu.aldex[order(abs(h.eu.aldex$effect),decreasing=TRUE),]
h.up.aldex <- h.up.aldex[order(abs(h.up.aldex$effect),decreasing=TRUE),]
h.extreme.aldex <- h.extreme.aldex[order(abs(h.extreme.aldex$effect),decreasing=TRUE),]

head(taxonomy[as.numeric(rownames(h.ep.aldex))])
# [1] "Bacteria;Bacteroidetes;Bacteroidia;Bacteroidales;Rikenellaceae;Alistipes;|100"               
# [2] "Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae;Coprococcus;|100"               
# [3] "Bacteria;Actinobacteria;Coriobacteriia;Coriobacteriales;Coriobacteriaceae;Adlercreutzia;|100"
# [4] "Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae;Incertae_Sedis;|97"             
# [5] "Bacteria;Firmicutes;Clostridia;Clostridiales;Ruminococcaceae;Subdoligranulum;|100"           
# [6] "Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae;Incertae_Sedis;|95"             
head(taxonomy[as.numeric(rownames(h.eu.aldex))])
# [1] "Bacteria;Bacteroidetes;Bacteroidia;Bacteroidales;Rikenellaceae;Alistipes;|100"                       
# [2] "Bacteria;Bacteroidetes;Bacteroidia;Bacteroidales;Rikenellaceae;Alistipes;|100"                       
# [3] "Bacteria;Firmicutes;Clostridia;Clostridiales;unclassified;unclassified;|94"                          
# [4] "Bacteria;Proteobacteria;Deltaproteobacteria;Desulfovibrionales;Desulfovibrionaceae;Bilophila;|100"   
# [5] "Bacteria;Firmicutes;Erysipelotrichia;Erysipelotrichales;Erysipelotrichaceae;Catenibacterium;|100"    
# [6] "Archaea;Euryarchaeota;Methanobacteria;Methanobacteriales;Methanobacteriaceae;Methanobrevibacter;|100"
head(taxonomy[as.numeric(rownames(h.up.aldex))])
# [1] "Bacteria;Firmicutes;Clostridia;Clostridiales;Ruminococcaceae;Faecalibacterium;|100"
# [2] "Bacteria;Firmicutes;Clostridia;Clostridiales;Ruminococcaceae;Faecalibacterium;|100"
# [3] "Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae;Incertae_Sedis;|97"   
# [4] "Bacteria;Bacteroidetes;Bacteroidia;Bacteroidales;Bacteroidaceae;Bacteroides;|100"  
# [5] "Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae;Coprococcus;|100"     
# [6] "Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae;Incertae_Sedis;|95"   
head(taxonomy[as.numeric(rownames(h.extreme.aldex))])
# [1] "Bacteria;Firmicutes;Clostridia;Clostridiales;Ruminococcaceae;Faecalibacterium;|100"        
# [2] "Bacteria;Firmicutes;Clostridia;Clostridiales;Ruminococcaceae;Faecalibacterium;|100"        
# [3] "Bacteria;Actinobacteria;Coriobacteriia;Coriobacteriales;Coriobacteriaceae;Collinsella;|100"
# [4] "Bacteria;Bacteroidetes;Bacteroidia;Bacteroidales;Bacteroidaceae;Bacteroides;|100"          
# [5] "Bacteria;Bacteroidetes;Bacteroidia;Bacteroidales;Prevotellaceae;Prevotella;|100"           
# [6] "Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae;Incertae_Sedis;|97" 

write.table(h.ep.aldex,file="explained_vs_protected_aldex_output.txt",sep="\t",quote=FALSE)
write.table(h.eu.aldex,file="explained_vs_unexplained_aldex_output.txt",sep="\t",quote=FALSE)
write.table(h.up.aldex,file="unexplained_vs_protected_aldex_output.txt",sep="\t",quote=FALSE)
write.table(h.extreme.aldex,file="top_decile_vs_bottoim_decile_aldex_output.txt",sep="\t",quote=FALSE)

# TODO: fix below

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
