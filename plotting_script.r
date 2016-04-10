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
names(otu.sp) <- colnames(otu.tab)

# get all samples with residual scores > 2 and < -2 (2 stdev away from median)
extreme <- otu.tab[which(metadata$Standardized.Residual >= 2 | metadata$Standardized.Residual <= -2),]
# get all samples with residual scores between 1 and 2 and -1 and -2
intermediate.top <- otu.tab[which(metadata$Standardized.Residual >= 1 & metadata$Standardized.Residual <2),]
intermediate.bottom <- otu.tab[which(metadata$Standardized.Residual <= -1 & metadata$Standardized.Residual > -2),]
intermediate <- intermediate <- t(data.frame(t(intermediate.top), t(intermediate.bottom), check.names=FALSE))
#separate by sex
extreme.m <- extreme[which(metadata[match(rownames(extreme),rownames(metadata)),"Sex"] == "M"),]
extreme.f <- extreme[which(metadata[match(rownames(extreme),rownames(metadata)),"Sex"] == "F"),]
intermediate.m <- intermediate[which(metadata[match(rownames(intermediate),rownames(metadata)),"Sex"] == "M"),]
intermediate.f <- intermediate[which(metadata[match(rownames(intermediate),rownames(metadata)),"Sex"] == "F"),]

make.groups <- function(mytable,metadata) {
  mygroups <- rep("NA",nrow(mytable))
  residuals <- metadata[match(rownames(mytable),rownames(metadata)),"Standardized.Residual"]
  mygroups[which(residuals <= -1)] <- "Protected"
  mygroups[which(residuals >= 1)] <- "Unexplained"
  mygroups[which(residuals < 1 & residuals > -1)] <- "Explained"
  return(as.factor(mygroups))
}

groups.extreme.m <- make.groups(extreme.m,metadata)
groups.extreme.f <- make.groups(extreme.f,metadata)
groups.intermediate.m <- make.groups(intermediate.m,metadata)
groups.intermediate.f <- make.groups(intermediate.f,metadata)

d.extreme.m <- t(extreme.m)
d.extreme.f <- t(extreme.f)
d.intermediate.m <- t(intermediate.m)
d.intermediate.f <- t(intermediate.f)

remove.zero.features <- function(mytable) {
  feature.sum <- apply(mytable,1,sum)
  mytable <- mytable[which(feature.sum > 0),]
  return(mytable)
}

d.extreme.m <- remove.zero.features(d.extreme.m)
d.extreme.f <- remove.zero.features(d.extreme.f)
d.intermediate.m <- remove.zero.features(d.intermediate.m)
d.intermediate.f <- remove.zero.features(d.intermediate.f)

# adjust zeros
d.extreme.m.adj.zero <- t(cmultRepl(t(d.extreme.m),method="CZM"))
d.extreme.f.adj.zero <- t(cmultRepl(t(d.extreme.f),method="CZM"))
d.intermediate.m.adj.zero <- t(cmultRepl(t(d.intermediate.m),method="CZM"))
d.intermediate.f.adj.zero <- t(cmultRepl(t(d.intermediate.f),method="CZM"))

d.extreme.m.names <- rownames(d.extreme.m.adj.zero)
d.extreme.f.names <- rownames(d.extreme.f.adj.zero)
d.intermediate.m.names <- rownames(d.intermediate.m.adj.zero)
d.intermediate.f.names <- rownames(d.intermediate.f.adj.zero)

taxa.col <- data.frame(otu.sp,otu.sp)
colnames(taxa.col) <- c("taxon","color")
taxa.col[,2] <- distinctColorPalette(length(taxa.col[,2]))

get.conds.df <- function(mygroups) {
  conds <- data.frame(as.character(mygroups))
  colnames(conds) <- "cond"
  return(conds)
}

conds.extreme.m <- get.conds.df(groups.extreme.m)
conds.extreme.f <- get.conds.df(groups.extreme.f)
conds.intermediate.m <- get.conds.df(groups.intermediate.m)
conds.intermediate.f <- get.conds.df(groups.intermediate.f)

palette=palette(c(rgb(1,0,0,0.6), "green","black", rgb(0,1,1,0.6)))

get.colors <- function(mygroups,dark="") {
  my.col <- rep("NA",length(mygroups))
  my.col[which(mygroups == "Explained")] <- paste(dark,"purple",sep="")
  my.col[which(mygroups == "Protected")] <- paste(dark,"blue",sep="")
  my.col[which(mygroups == "Unexplained")] <- paste(dark,"red",sep="")
  return(my.col)
}

col.extreme.m <- get.colors(groups.extreme.m,dark="dark")
col.extreme.f <- get.colors(groups.extreme.f,dark="dark")
col.intermediate.m <- get.colors(groups.intermediate.m)
col.intermediate.f <- get.colors(groups.intermediate.f)

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

d.extreme.m.clr <- plot.biplot(d.extreme.m.adj.zero,col.extreme.m,groups.extreme.m,c("darkblue","darkred"))
d.extreme.f.clr <- plot.biplot(d.extreme.f.adj.zero,col.extreme.f,groups.extreme.f,c("darkblue","darkred"))
d.intermediate.m.clr <- plot.biplot(d.intermediate.m.adj.zero,col.intermediate.m,groups.intermediate.m,c("blue","red"))
d.intermediate.f.clr <- plot.biplot(d.intermediate.f.adj.zero,col.intermediate.f,groups.intermediate.f,c("blue","red"))

dev.off()

plot.dendogram.barplot <- function(d.clr, d.adj.zero, taxa.col, d.names, groups) {
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
  barplot(d.acomp, legend.text=F, col=as.character(taxa.col[match(colnames(d.acomp),rownames(taxa.col)),2]), axisnames=F, border=NA, xpd=T)
  par(mar=c(0,1,1,1)+0.1)
  # and the legend
  plot(1,2, pch = 1, lty = 1, ylim=c(-20,20), type = "n", axes = FALSE, ann = FALSE)
  legend(x="center", legend=taxa.col[match(d.names,rownames(taxa.col)),1], col=as.character(taxa.col[match(d.names,rownames(taxa.col)),2]), lwd=5, cex=.3, border=NULL,ncol=3)
}

pdf("dendogram_barplot.pdf")

plot.dendogram.barplot(d.extreme.m.clr, d.extreme.m.adj.zero, taxa.col,d.extreme.m.names,groups.extreme.m)
plot.dendogram.barplot(d.extreme.f.clr, d.extreme.f.adj.zero, taxa.col,d.extreme.f.names,groups.extreme.f)
plot.dendogram.barplot(d.intermediate.m.clr, d.intermediate.m.adj.zero, taxa.col,d.intermediate.m.names,groups.intermediate.m)
plot.dendogram.barplot(d.intermediate.f.clr, d.intermediate.f.adj.zero, taxa.col,d.intermediate.f.names,groups.intermediate.f)

d.clr <- d.extreme.f.clr
d.names <- d.extreme.f.names
groups <- groups.extreme.f
d.adj.zero <- d.extreme.f.adj.zero


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
# [1] -0.385572
plot(h.eu.aldex$effect, h.up.aldex$effect, pch=19,col=mycolor, main="Effect sizes of Explained vs. Unexplained\nand Unexplained vs. Protected",xlab="Explained vs. Unexplained",ylab="Unexplained vs. Protected")
cor(h.eu.aldex$effect, y = h.up.aldex$effect, use = "everything", method = "spearman")
# [1] 0.6208453
plot(h.up.aldex$effect, h.extreme.aldex$effect, pch=19,col=mycolor, main="Effect sizes of Unexplained vs. Protected\nand Top decile vs. Bottom decile subset",xlab="Unexplained vs. Protected",ylab="Top decile vs. Bottom decile subset")
cor(h.up.aldex$effect, y = h.extreme.aldex$effect, use = "everything", method = "spearman")
# [1] 0.5096007
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
write.table(h.extreme.aldex,file="top_decile_vs_bottom_decile_aldex_output.txt",sep="\t",quote=FALSE)

# make plots with top and bottom decile effect sizes from extreme comparison colored
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
h.extreme.aldex <- h.extreme.aldex[order(h.extreme.aldex$effect,decreasing=TRUE),]
h.ep.aldex <- h.ep.aldex[match(rownames(h.extreme.aldex), rownames(h.ep.aldex)),]
h.eu.aldex <- h.eu.aldex[match(rownames(h.extreme.aldex), rownames(h.eu.aldex)),]
h.up.aldex <- h.up.aldex[match(rownames(h.extreme.aldex), rownames(h.up.aldex)),]

effect.groups <- rep("NA",nrow(h.extreme.aldex))
effect.groups[c(1:decile)] <- "Top decile"
effect.groups[c((decile+1):(length(effect.groups)-decile))] <- "Middle half"
effect.groups[c((length(effect.groups)-decile+1):length(effect.groups))] <- "Bottom decile"
effect.groups <- as.factor(effect.groups)

plot(h.ep.aldex$effect, h.up.aldex$effect, pch=19,col=effect.groups, main="Effect sizes of explained vs protected\nand unexplained vs protected",xlab="explained vs protected",ylab="unexplained vs protected")
cor(h.ep.aldex$effect, y = h.up.aldex$effect, use = "everything", method = "spearman")
# [1] -0.385572

plot(h.eu.aldex$effect, h.up.aldex$effect, pch=19,col=effect.groups, main="Effect sizes of explained vs unexplained\nand unexplained vs protected",xlab="explained vs unexplained",ylab="unexplained vs protected")
cor(h.eu.aldex$effect, y = h.up.aldex$effect, use = "everything", method = "spearman")
# [1] 0.6208453

plot(h.up.aldex$effect, h.extreme.aldex$effect, pch=19,col=effect.groups, main="Effect sizes of unexplained vs protected\nand extreme deciles",xlab="unexplained vs protected",ylab="extreme deciles")
cor(h.up.aldex$effect, y = h.extreme.aldex$effect, use = "everything", method = "spearman")
# [1] 0.5096007