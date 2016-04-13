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

pdf("aldex_plots.pdf")

plot.aldex.volcano(data.frame(d.extreme.m),as.character(groups.extreme.m))
plot.aldex.volcano(data.frame(d.extreme.f),as.character(groups.extreme.f))
plot.aldex.volcano(data.frame(d.intermediate.m),as.character(groups.intermediate.m))
plot.aldex.volcano(data.frame(d.intermediate.f),as.character(groups.intermediate.f))

dev.off()

extreme.m.aldex <- aldex(data.frame(d.extreme.m),as.character(groups.extreme.m))
extreme.f.aldex <- aldex(data.frame(d.extreme.f),as.character(groups.extreme.f))
intermediate.m.aldex <- aldex(data.frame(d.intermediate.m),as.character(groups.intermediate.m))
intermediate.f.aldex <- aldex(data.frame(d.intermediate.f),as.character(groups.intermediate.f))

mycolor <- c(col2rgb("turquoise4"))
red <- mycolor[1]
green <- mycolor[2]
blue <- mycolor[3]
mycolor <- rgb(red/255, green/255, blue/255, 0.3)

pdf("16S_effect_sizes.pdf")
common.features <- intersect(intersect(intersect(rownames(extreme.m.aldex),rownames(extreme.f.aldex)), rownames(intermediate.m.aldex)), rownames(intermediate.f.aldex))
plot(extreme.m.aldex[match(common.features,rownames(extreme.m.aldex)),"effect"], extreme.f.aldex[match(common.features,rownames(extreme.f.aldex)),"effect"], pch=19,col=mycolor, main="Effect sizes of OTUs between male\nand female extreme residual scorers",xlab="Male sex",ylab="Female sex")
cor(extreme.m.aldex[match(common.features,rownames(extreme.m.aldex)),"effect"], y = extreme.f.aldex[match(common.features,rownames(extreme.f.aldex)),"effect"], use = "everything", method = "spearman")
# [1] -0.08417705
plot(intermediate.m.aldex[match(common.features,rownames(intermediate.m.aldex)),"effect"], intermediate.f.aldex[match(common.features,rownames(intermediate.f.aldex)),"effect"], pch=19,col=mycolor, main="Effect sizes of OTUs between male\nand female intermediate residual scorers",xlab="Male sex",ylab="Female sex")
cor(intermediate.m.aldex[match(common.features,rownames(intermediate.m.aldex)),"effect"], y = intermediate.f.aldex[match(common.features,rownames(intermediate.f.aldex)),"effect"], use = "everything", method = "spearman")
# [1] 0.2434751
plot(extreme.m.aldex[match(common.features,rownames(extreme.m.aldex)),"effect"], intermediate.m.aldex[match(common.features,rownames(intermediate.m.aldex)),"effect"], pch=19,col=mycolor, main="Effect sizes of OTUs between male\nextreme and intermediate residual scorers",xlab="Extreme residual scores",ylab="Intermediate residual scores")
cor(extreme.m.aldex[match(common.features,rownames(extreme.m.aldex)),"effect"], y = intermediate.m.aldex[match(common.features,rownames(intermediate.m.aldex)),"effect"], use = "everything", method = "spearman")
# [1] 0.008379229
plot(extreme.f.aldex[match(common.features,rownames(extreme.f.aldex)),"effect"], intermediate.f.aldex[match(common.features,rownames(intermediate.f.aldex)),"effect"], pch=19,col=mycolor, main="Effect sizes of OTUs between female\nextreme and intermediate residual scorers",xlab="Extreme residual scores",ylab="Intermediate residual scores")
cor(extreme.f.aldex[match(common.features,rownames(extreme.f.aldex)),"effect"], y = intermediate.f.aldex[match(common.features,rownames(intermediate.f.aldex)),"effect"], use = "everything", method = "spearman")
# [1] -0.05444832
dev.off()

# MALE VS FEMALE ALDEX PLOTS

d <- data.frame(matrix(nrow=nrow(taxa.col), ncol=(ncol(d.extreme.m) + ncol(d.extreme.f) + ncol(d.intermediate.m) + ncol(d.intermediate.f))))
rownames(d) <- rownames(taxa.col)
d[match(rownames(d.extreme.m),rownames(d)),c(1:ncol(d.extreme.m))] <- d.extreme.m
colnames(d)[c(1:ncol(d.extreme.m))] <- colnames(d.extreme.m)
index <- ncol(d.extreme.m) + 1
d[match(rownames(d.extreme.f),rownames(d)),c(index:(index + ncol(d.extreme.f) - 1))] <- d.extreme.f
colnames(d)[c(index:(index + ncol(d.extreme.f) - 1))] <- colnames(d.extreme.f)
index <- index + ncol(d.extreme.f)
d[match(rownames(d.intermediate.m),rownames(d)),c(index:(index + ncol(d.intermediate.m) - 1))] <- d.intermediate.m
colnames(d)[c(index:(index + ncol(d.intermediate.m) - 1))] <- colnames(d.intermediate.m)
index <- index + ncol(d.intermediate.m)
d[match(rownames(d.intermediate.f),rownames(d)),c(index:(index + ncol(d.intermediate.f) - 1))] <- d.intermediate.f
colnames(d)[c(index:(index + ncol(d.intermediate.f) - 1))] <- colnames(d.intermediate.f)

d.conds <- c(rep("Male",ncol(d.extreme.m)),rep("Female",ncol(d.extreme.f)),rep("Male",ncol(d.intermediate.m)),rep("Female",ncol(d.intermediate.f)))

d.aldex <- aldex(d,d.conds)


pdf("aldex_male_vs_female.pdf")

aldex.plot(d.aldex,type="MA")
aldex.plot(d.aldex,type="MW")

dev.off()

d.aldex <- d.aldex[order(abs(d.aldex$effect),decreasing=TRUE),]
write.table(d.aldex,file="ALDEx_male_vs_female_output.txt",sep="\t",quote=FALSE)

unique.taxa <- paste(rownames(taxa.col),taxa.col[,1],sep=":")
# the cutoff of 0.4 is chosen so that there is at least one qualifying taxa in each group
print("Unique taxa with effect size >= 0.4 in Extreme M condition:")
print(paste(unique.taxa[match(rownames(extreme.m.aldex)[which(extreme.m.aldex$effect >= 0.4)],rownames(taxa.col))],collapse="\n"))
print("Unique taxa with effect size >= 0.4 in Extreme F condition:")
print(paste(unique.taxa[match(rownames(extreme.f.aldex)[which(extreme.f.aldex$effect >= 0.4)],rownames(taxa.col))],collapse="\n"))
print("Unique taxa with effect size >= 0.4 in Intermediate M condition:")
print(paste(unique.taxa[match(rownames(intermediate.m.aldex)[which(intermediate.m.aldex$effect >= 0.4)],rownames(taxa.col))],collapse="\n"))
print("Unique taxa with effect size >= 0.4 in Intermediate F condition:")
print(paste(unique.taxa[match(rownames(intermediate.f.aldex)[which(intermediate.f.aldex$effect >= 0.4)],rownames(taxa.col))],collapse="\n"))

# 4 common elements in "male, extreme" and "female, extreme":
# 52:Adlercreutzia|100
# 94:Coprococcus|100
# 101:Eggerthella|100
# 734:Incertae_Sedis|90

print("Unique taxa with effect size <= -0.4 in Extreme M condition:")
print(paste(unique.taxa[match(rownames(extreme.m.aldex)[which(extreme.m.aldex$effect <= -0.4)],rownames(taxa.col))],collapse="\n"))
print("Unique taxa with effect size <= -0.4 in Extreme F condition:")
print(paste(unique.taxa[match(rownames(extreme.f.aldex)[which(extreme.f.aldex$effect <= -0.4)],rownames(taxa.col))],collapse="\n"))
print("Unique taxa with effect size <= -0.4 in Intermediate M condition:")
print(paste(unique.taxa[match(rownames(intermediate.m.aldex)[which(intermediate.m.aldex$effect <= -0.4)],rownames(taxa.col))],collapse="\n"))
print("Unique taxa with effect size <= -0.4 in Intermediate F condition:")
print(paste(unique.taxa[match(rownames(intermediate.f.aldex)[which(intermediate.f.aldex$effect <= -0.4)],rownames(taxa.col))],collapse="\n"))

# 2 common elements in "male, extreme" and "female, extreme":
# 85:Incertae_Sedis|84
# 1855:Incertae_Sedis|99

# 4 common elements in "male, extreme" and "female, intermed":
# 96:Faecalibacterium|100
# 44:Roseburia|94
# 1223:Incertae_Sedis|93
# 783:Ruminococcus|100

# 2 common elements in "male, extreme", "female, extreme" and "female, intermed":
# 791:Faecalibacterium|100
# 645:Incertae_Sedis|97

# 1 common element in "female, extreme" and "female, intermed":
# 68:Incertae_Sedis|97

# intermediate effects dont seem to tend to be higher than extreme effects
summary(extreme.m.aldex$effect)
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -0.98160 -0.19640 -0.01633 -0.02153  0.16320  0.76090 
summary(extreme.f.aldex$effect)
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -0.75440 -0.14740  0.02674  0.02125  0.18290  0.98720 
summary(intermediate.m.aldex$effect)
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -0.51140 -0.12150 -0.01007 -0.02019  0.09708  0.45800 
summary(intermediate.f.aldex$effect)
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# -0.684500 -0.127700 -0.003546 -0.010050  0.110800  0.559700 

summary(extreme.m.aldex[match(common.features,rownames(extreme.m.aldex)),"effect"] - intermediate.m.aldex[match(common.features,rownames(intermediate.m.aldex)),"effect"])
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# -0.846600 -0.228200 -0.016500 -0.000976  0.205000  0.835600 
summary(extreme.f.aldex[match(common.features,rownames(extreme.f.aldex)),"effect"] - intermediate.f.aldex[match(common.features,rownames(intermediate.f.aldex)),"effect"])
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -0.81830 -0.23020  0.01923  0.03125  0.26710  1.02800 


# EXPLORE PREVOTELLA
# One Prevotella is missing in extreme.f group. See what it's up to in the other groups.
prevotella <- rownames(taxa.col)[which(!(rownames(taxa.col) %in% rownames(d.extreme.f)))]
print(paste("Effect size of missing prevotella in Extreme M:",extreme.m.aldex$effect[which(rownames(extreme.m.aldex) ==  prevotella)]))
# [1] "Effect size of missing prevotella in Extreme M: -0.112519887389879"
print(paste("Effect size of missing prevotella in Intermediate M:",intermediate.m.aldex$effect[which(rownames(intermediate.m.aldex) ==  prevotella)]))
# [1] "Effect size of missing prevotella in Intermediate M: -0.00244435674671425"
print(paste("Effect size of missing prevotella in Intermediate F:",intermediate.f.aldex$effect[which(rownames(intermediate.f.aldex) ==  prevotella)]))
# [1] "Effect size of missing prevotella in Intermediate F: -0.0249854285733147"

# SHANNON DIVERSITY
library(vegan)
diversity.extreme.m <- diversity(t(d.extreme.m),index="shannon")
diversity.extreme.f <- diversity(t(d.extreme.f),index="shannon")
diversity.intermediate.m <- diversity(t(d.intermediate.m),index="shannon")
diversity.intermediate.f <- diversity(t(d.intermediate.f),index="shannon")

# protected vs unexplained
print(paste("Effect size of protected vs. unexplained diversity for male sex with extreme residuals:", (abs(mean(diversity.extreme.m[which(groups.extreme.m == "Protected")]) - mean(diversity.extreme.m[which(groups.extreme.m == "Unexplained")])) / sd(diversity.extreme.m))))
# [1] "Effect size of protected vs. unexplained diversity for male sex with extreme residuals: 0.244250455280845"
print(paste("Effect size of protected vs. unexplained diversity for female sex with extreme residuals:", (abs(mean(diversity.extreme.f[which(groups.extreme.f == "Protected")]) - mean(diversity.extreme.f[which(groups.extreme.f == "Unexplained")])) / sd(diversity.extreme.f))))
# [1] "Effect size of protected vs. unexplained diversity for female sex with extreme residuals: 0.222511433991768"
print(paste("Effect size of protected vs. unexplained diversity for male sex with intermediate residuals:", (abs(mean(diversity.intermediate.m[which(groups.intermediate.m == "Protected")]) - mean(diversity.intermediate.m[which(groups.intermediate.m == "Unexplained")])) / sd(diversity.intermediate.m))))
# [1] "Effect size of protected vs. unexplained diversity for male sex with intermediate residuals: 0.276879902237068"
print(paste("Effect size of protected vs. unexplained diversity for female sex with intermediate residuals:", (abs(mean(diversity.intermediate.f[which(groups.intermediate.f == "Protected")]) - mean(diversity.intermediate.f[which(groups.intermediate.f == "Unexplained")])) / sd(diversity.intermediate.f))))
# [1] "Effect size of protected vs. unexplained diversity for female sex with intermediate residuals: 0.2381063694763"

# extreme vs intermediate
print(paste("Effect size of extreme vs. intermediate diversity for male sex:", (abs(mean(diversity.extreme.m) - mean(diversity.intermediate.m)) / sd(c(diversity.extreme.m, diversity.intermediate.m)))))
# [1] "Effect size of extreme vs. intermediate diversity for male sex: 0.556190597796183"
print(paste("Effect size of extreme vs. intermediate diversity for female sex:", (abs(mean(diversity.extreme.f) - mean(diversity.intermediate.f)) / sd(c(diversity.extreme.f, diversity.intermediate.f)))))
# [1] "Effect size of extreme vs. intermediate diversity for female sex: 0.238379974990715"

# male vs female
print(paste("Effect size of male vs. female sex for extreme residuals:", (abs(mean(diversity.extreme.m) - mean(diversity.extreme.f)) / sd(c(diversity.extreme.m, diversity.extreme.f)))))
# [1] "Effect size of male vs. female sex for extreme residuals: 0.297710846323797"
print(paste("Effect size of male vs. female sex for intermediate residuals:", (abs(mean(diversity.intermediate.m) - mean(diversity.intermediate.f)) / sd(c(diversity.intermediate.m, diversity.intermediate.f)))))
# [1] "Effect size of male vs. female sex for intermediate residuals: 0.00689408376308769"

plot.extremes.metadata.aldex <- function(otu.tab,metadata.category,metadataname) {
	residuals <- order(metadata.category,decreasing=TRUE)
	decile <- round(length(residuals)/10)
	top <- residuals[c(1:decile)]
	bottom <- residuals[c((length(residuals)-decile + 1):length(residuals))]
	otu.tab.metadata <- otu.tab[c(top,bottom),]
	groups <- c(rep(paste("Top",metadataname),length(top)),rep(paste("Bottom",metadataname),length(bottom)))
  
  plot.aldex.volcano(data.frame(t(otu.tab.metadata)),groups)
}

pdf("aldex_plots_metadata_category_extremes.pdf")
plot.extremes.metadata.aldex(otu.tab,metadata$TPAmm2,"TPA mm2")
plot.extremes.metadata.aldex(otu.tab,metadata$Age,"Age")
plot.extremes.metadata.aldex(otu.tab,metadata$PackYears,"Pack Years")
plot.extremes.metadata.aldex(otu.tab,metadata$BPSys,"Systolic BP")
plot.extremes.metadata.aldex(otu.tab,metadata$BPDias,"Diastolic BP")
plot.extremes.metadata.aldex(otu.tab,metadata$chol,"cholesteral")
plot.extremes.metadata.aldex(otu.tab,metadata$trig,"triglyceride")
plot.extremes.metadata.aldex(otu.tab,metadata$hdl,"HDL")
plot.extremes.metadata.aldex(otu.tab,metadata$ldl,"LDL")
plot.extremes.metadata.aldex(otu.tab,metadata$Standardized.Predicted.Value,"Predicted")
dev.off()

#standardized predicted value has some significants - explore
spv <- order(metadata$Standardized.Predicted.Value,decreasing=TRUE)
decile <- round(length(spv)/10)
top <- spv[c(1:decile)]
bottom <- spv[c((length(spv)-decile + 1):length(spv))]
otu.tab.svp <- otu.tab[c(top,bottom),]
groups <- c(rep(paste("Top","Predicted"),length(top)),rep(paste("Bottom","Predicted"),length(bottom)))

spv.aldex <- aldex(data.frame(t(otu.tab.svp)),groups)
spv.aldex <- spv.aldex[order(abs(spv.aldex$effect),decreasing=TRUE),]
spv.otu <- rownames(spv.aldex)[which(spv.aldex$we.eBH <= 0.05)]
spv.taxa <- taxa.col[match(spv.otu,rownames(taxa.col)),1]
print(spv.taxa)
# [1] unclassified|100 Coprococcus|91  
# 138 Levels: Acidaminococcus|100 Adlercreutzia|100 ... Weissella|100
print(extreme.m.aldex[match(spv.otu,rownames(extreme.m.aldex)),])
# rab.all rab.win.Protected rab.win.Unexplained  diff.btw diff.win
# 86 0.8902948         0.6899341            1.259840 0.7372637 2.727213
# 33 4.6199940         4.5259409            4.839426 0.1399725 1.820813
# diff.btw.025 diff.btw.975 diff.win.025 diff.win.975     effect effect.025
# 86    -3.890029     6.002033    0.3194901     5.541354 0.27594571  -4.119158
# 33    -3.473844     2.987637    0.2815473     3.857473 0.07190766  -4.237978
# effect.975   overlap     we.ep    we.eBH     wi.ep    wi.eBH
# 86   7.189024 0.3281251 0.4185800 0.8736268 0.2348797 0.7999611
# 33   3.815769 0.4687500 0.9171129 0.9879872 0.8926024 0.9865305
print(extreme.f.aldex[match(spv.otu,rownames(extreme.f.aldex)),])
# rab.all rab.win.Protected rab.win.Unexplained    diff.btw diff.win
# 86 0.8453854          1.283060           0.5652007 -0.62050192 2.186366
# 33 4.4924272          4.490632           4.4949703 -0.02118608 1.917428
# diff.btw.025 diff.btw.975 diff.win.025 diff.win.975       effect effect.025
# 86    -3.868434     3.273848    0.2707417     4.758016 -0.357925670  -5.655260
# 33    -4.801129     3.223706    0.2533891     6.020941 -0.007999378  -5.426252
# effect.975   overlap     we.ep    we.eBH     wi.ep    wi.eBH
# 86   3.856239 0.3318486 0.3368028 0.9227263 0.3501991 0.9326479
# 33   4.017823 0.4922049 0.8432237 0.9905316 0.9258278 0.9922410
print(intermediate.m.aldex[match(spv.otu,rownames(intermediate.m.aldex)),])
# rab.all rab.win.Protected rab.win.Unexplained    diff.btw diff.win
# 86 1.346272          1.441525            1.267438 -0.03775829 2.696193
# 33 4.912891          4.866859            4.955281  0.01898387 1.776607
# diff.btw.025 diff.btw.975 diff.win.025 diff.win.975       effect effect.025
# 86    -4.353421     4.722191    0.4328533     5.589769 -0.013551302  -3.328893
# 33    -3.267390     3.263029    0.3085356     3.955361  0.008311554  -3.683919
# effect.975   overlap     we.ep    we.eBH     wi.ep    wi.eBH
# 86   3.635860 0.4925272 0.8729045 0.9739340 0.8055559 0.9693421
# 33   3.337985 0.4952446 0.9292273 0.9846214 0.9162946 0.9874256
print(intermediate.f.aldex[match(spv.otu,rownames(intermediate.f.aldex)),])
# rab.all rab.win.Protected rab.win.Unexplained    diff.btw diff.win
# 86 1.396652          1.297878            1.438710 -0.03533186 2.408819
# 33 4.123899          4.111840            4.151104 -0.08414671 1.554466
# diff.btw.025 diff.btw.975 diff.win.025 diff.win.975      effect effect.025
# 86    -4.673503     4.404538    0.3646155     5.140760 -0.01577379  -4.668171
# 33    -2.818697     2.718209    0.2414507     3.738445 -0.05109171  -3.654728
# effect.975   overlap     we.ep    we.eBH     wi.ep    wi.eBH
# 86   3.873488 0.4903581 0.8382479 0.9780627 0.8288777 0.9656435
# 33   3.940651 0.4756658 0.9055308 0.9850980 0.8700055 0.9748175

# STENOSIS RESIDUALS
stenosis.metadata <- read.table("data/stenosis_residuals.tsv",sep="\t",check.names=FALSE,quote="",header=TRUE,row.names=1)
colnames(stenosis.metadata) <- gsub(" ",".",colnames(stenosis.metadata))
stenosis.svp <- stenosis.metadata$Standardized.Predicted.Value[match(metadata$idno,stenosis.metadata$idno)]
stenosis.res <- stenosis.metadata$Standardized.Residual[match(metadata$idno,stenosis.metadata$idno)]
pdf("aldex_plots_stenosis.pdf")
plot.extremes.metadata.aldex(otu.tab,stenosis.svp,"Predicted")
plot.extremes.metadata.aldex(otu.tab,stenosis.res,"Residual")
plot(metadata$Standardized.Residual,stenosis.res)
dev.off()

# STENOSIS EFFECT SIZES VS TPA EFFECT SIZES
stenosis.order <- order(stenosis.res,decreasing=TRUE)
residual.order <- order(metadata$Standardized.Residual,decreasing=TRUE)
decile <- round(length(stenosis.order)/10)

stenosis.top <- rownames(metadata)[stenosis.order[c(1:decile)]]
residual.top <- rownames(metadata)[residual.order[c(1:decile)]]
stenosis.bottom <- rownames(metadata)[stenosis.order[c((length(stenosis.order) - decile + 1):length(stenosis.order))]]
residual.bottom <- rownames(metadata)[residual.order[c((length(residual.order) - decile + 1):length(residual.order))]]

intersecting.samples <- c(intersect(stenosis.top, residual.top), intersect(stenosis.bottom, residual.bottom), intersect(stenosis.top, residual.bottom), intersect(stenosis.bottom, residual.top))
stenosis.top <- stenosis.top[which(!(stenosis.top %in% intersecting.samples))]
stenosis.bottom <- stenosis.bottom[which(!(stenosis.bottom %in% intersecting.samples))]
residual.top <- residual.top[which(!(residual.top %in% intersecting.samples))]
residual.bottom <- residual.bottom[which(!(residual.bottom %in% intersecting.samples))]

stenosis.counts <- data.frame(otu.tab[c(match(stenosis.top,rownames(otu.tab)),match(stenosis.bottom,rownames(otu.tab))),],check.names=FALSE)
residual.counts <- data.frame(otu.tab[c(match(residual.top,rownames(otu.tab)),match(residual.bottom,rownames(otu.tab))),],check.names=FALSE)
stenosis.conds <- c(rep("Top decile",length(stenosis.top)),rep("Bottom decile",length(stenosis.bottom)))
residual.conds <- c(rep("Top decile",length(residual.top)),rep("Bottom decile",length(residual.bottom)))
stenosis.aldex <- aldex(data.frame(t(stenosis.counts)),stenosis.conds)
residual.aldex <- aldex(data.frame(t(residual.counts)),residual.conds)

mycolor <- c(col2rgb("turquoise4"))
red <- mycolor[1]
green <- mycolor[2]
blue <- mycolor[3]
mycolor <- rgb(red/255, green/255, blue/255, 0.3)

pdf("stenosis_effect_sizes_comparison.pdf")

plot(residual.aldex$effect, stenosis.aldex$effect, pch=19,col=mycolor, main="Effect sizes of OTUs between stenosis\nand TPA mm2 extreme decile residual scorers",xlab="TPA mm2",ylab="Stenosis residuals")
cor(residual.aldex$effect, y = stenosis.aldex$effect, use = "everything", method = "spearman")
# [1] -0.06968418

dev.off()
