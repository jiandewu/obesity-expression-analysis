library(limma)
library(edgeR)


options(stringsAsFactors = FALSE)

data_dir <- Sys.getenv("OBESITY_DATA_DIR", unset = "data")
rawTNBC <- read.table(file.path(data_dir, "GSE76124", "03_data_exp_genes_165_bmisorted.txt"), header = TRUE, sep="\t")
rawTNBC[6,] <- "TNBC"
rawTNBC <- rawTNBC[c(-2,-5),-2]
rawTNBC[1,1] <- "!Characteristics_Age"
rawTNBC[2,1] <- "!Characteristics_Race"
rawTNBC[3,1] <- "!Characteristics_BMI"
rawTNBC[4,1] <- "!Group"
rawTNBC[1,] <- gsub("age \\(years\\): ", "", rawTNBC[1,])
rawTNBC[2,] <- gsub("race: ", "", rawTNBC[2,])
rawTNBC[2,] <- gsub("Caucasian", "white", rawTNBC[2,])
rawTNBC[2,] <- gsub("Asian/Pacific Islander", "Asian", rawTNBC[2,])
# rawTNBC$GSM1974705[2] <- "NA"
rawTNBC[3,] <- gsub("body mass index: ", "", rawTNBC[3,])

rawCtrl <- read.table(file.path(data_dir, "GSE102088", "03_filtered.csv"), header = TRUE, sep=",")
rawCtrl[4,] <- "Normal"
rawCtrl <- rawCtrl[,-1]
rawCtrl[1,1] <- "!Characteristics_Age"
rawCtrl[2,1] <- "!Characteristics_Race"
rawCtrl[3,1] <- "!Characteristics_BMI"
rawCtrl[4,1] <- "!Group"
rawCtrl[1,] <- gsub("age: ", "", rawCtrl[1,])
rawCtrl[2,] <- gsub("race: ", "", rawCtrl[2,])
rawCtrl[3,] <- gsub("bmi: ", "", rawCtrl[3,])
rawCtrl[,1] <- gsub(" ", "", rawCtrl[,1])
colnames(rawCtrl)[1] <- "X1"

rawdata <- merge(rawTNBC,rawCtrl,by = "X1")

rawdata[c(-1:-4),-1] <- round(2^(sapply(rawdata[c(-1:-4),-1], as.numeric)))

# Rename duplicate rows
row.names(rawdata) <- ave(as.character(rawdata[,1]), rawdata[,1], FUN=function(x) if (length(x)>1) paste0(x[1], '(', seq_along(x), ')') else x[1])

data1 <- DGEList(counts = sapply(rawdata[c(-1:-4),-1], as.numeric))

bmi_numeric <- as.numeric(unlist(rawdata[2,-1], use.names = FALSE))
bmi <- cut(
  bmi_numeric,
  breaks = c(-Inf, 24.99, 29.99, Inf),
  labels = c("Lean", "OverWeight", "Obese")
)
racelevel <- rawdata[3,-1]
race <- factor(unlist(racelevel, use.names = FALSE))
grouplevel <- rawdata[4,-1]
group <- factor(unlist(grouplevel, use.names = FALSE))
data1$samples$bmi <- bmi
data1$samples$race <- race
data1$samples$group <- group

table(data1$samples$group,data1$samples$bmi,data1$samples$race)

# Transformations from the raw-scale
seqdepth_mean <- mean(data1$samples$lib.size) * 1e-6
seqdepth_median <- median(data1$samples$lib.size) * 1e-6

summary(data1$samples$lib.size * 1e-6)
nsamples <- ncol(data1)

table(rowSums(data1$counts==0)==nsamples)

# Plotting the distribution log-CPM values shows that
# a sizeable proportion of genes within each sample are
# either unexpressed or lowly-expressed with log-CPM values
# that are small or negative

lcpm.cutoff <- log2(10/seqdepth_median + 2/seqdepth_mean)
library(RColorBrewer)
col <- colorRampPalette(brewer.pal(12, "Paired"))(nsamples)

lcpm <- cpm(data1, log=TRUE, prior.count=2)
par(mfrow=c(1,2))
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.4), las=2, main="", xlab="")
title(main="A. Raw data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", colnames(data1), text.col=col, bty="n")

# Removing genes that are lowly expressed
keep.exprs <- filterByExpr(data1, group=group)
data2 <- data1[keep.exprs,, keep.lib.sizes=FALSE]
dim(data2)
## [1] 792 279

lcpm <- cpm(data2, log=TRUE, prior.count=2)
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.4), las=2, main="", xlab="")
title(main="B. Filtered data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", colnames(data1), text.col=col, bty="n")

# Normalising gene expression distributions
par(mfrow=c(1,2))
lcpm <- cpm(data2, log=TRUE, prior.count=2)
boxplot(lcpm, las=2, col=col, main="")
title(main="A. Unnormalised data",ylab="Log-cpm")
data2 <- calcNormFactors(data2, method = "TMM")
data2$samples$norm.factors
lcpm <- cpm(data2, log=TRUE, prior.count=2)
boxplot(lcpm, las=2, col=col, main="")
title(main="B. Normalised data",ylab="Log-cpm")

# Unsupervised clustering of samples
lcpm <- cpm(data2, log=TRUE, prior.count=2)
par(mfrow=c(1,2))
col.group <- group
levels(col.group) <- brewer.pal(max(3, nlevels(col.group)), "Set1")[seq_len(nlevels(col.group))]
col.group <- as.character(col.group)
col.bmi <- bmi
levels(col.bmi) <- brewer.pal(max(3, nlevels(col.bmi)), "Set2")[seq_len(nlevels(col.bmi))]
col.bmi <- as.character(col.bmi)
plotMDS(lcpm, labels=group, col=col.group)
title(main="A. BC Category")
plotMDS(lcpm, labels=bmi, col=col.bmi, dim=c(3,4))
title(main="B. BMI")


design <- model.matrix(~0 + bmi)
colnames(design) <- sub("^bmi", "", colnames(design))
contr.matrix <- makeContrasts(
  Lean.vs.Obese = Lean - Obese,
  Lean.vs.OverWeight = Lean - OverWeight,
  levels = colnames(design)
)
contr.matrix


par(mfrow=c(1,2))
v <- voom(data2, design, plot=TRUE)
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
plotSA(efit, main="Final model: Mean-variance trend")
dt <- decideTests(efit)
summary(dt)

# tfit <- treat(vfit, lfc=1)
# summary(decideTests(tfit))


plotMD(efit, column=1, status=dt[,1], main=colnames(efit)[1], xlim=c(-3,13))

# vennDiagram(dt[,c(1,4)], circle.col=c("turquoise", "salmon"))
#
# write.fit(efit, dt, file="results.txt")

# de.common <- which(dt[,1]!=0 & dt[,4]!=0)
de.Lean.vs.Obese <- which(dt[,1] != 0)
Lean.vs.Obese <- topTable(efit, coef=1, n=Inf)
up.gene <- row.names(Lean.vs.Obese[Lean.vs.Obese$adj.P.Val<0.05&Lean.vs.Obese$logFC>0,])
down.gene <- row.names(Lean.vs.Obese[Lean.vs.Obese$adj.P.Val<0.05&Lean.vs.Obese$logFC<0,])

# Volcano draw
# volcanoplot(efit, coef=4,  main = "Lean.vs.Obese", highlight=8, names = row.names(efit$p.value), hl.col="blue", xlim=c(-3,3), ylim=c(0,5))
# abline(v=c(-1,1), h=3, lty=3)
library(ggplot2)
logFC.cutoff <- 1
adj.P.Val.cutoff <- 0.05
Volcano.color <- as.factor(ifelse(Lean.vs.Obese$adj.P.Val<adj.P.Val.cutoff &
                                    abs(Lean.vs.Obese$logFC)>=logFC.cutoff,
                                  ifelse(Lean.vs.Obese$logFC > logFC.cutoff,'Up','Down'),'Not'))

p <- ggplot(data=Lean.vs.Obese, aes(x=logFC, y=-log10(adj.P.Val), col=Volcano.color)) +
  geom_point(alpha=0.8, size = 0.8) +
  theme_bw(base_size = 15) +
  labs(x="log2 (fold change)",y="-log10 (adj.p-value)",title = "Lean.vs.Obese") +
  geom_hline(yintercept = -log10(adj.P.Val.cutoff), lty=4,col="black",lwd=0.6) +
  geom_vline(xintercept = c(-logFC.cutoff, logFC.cutoff), lty=4,col="black",lwd=0.6) +
  scale_color_manual(name = "", values = c("green", "red", "gray"), limits = c("Up", "Down", "Not")) +
  theme(panel.grid = element_blank(),
        plot.title = element_text(color="black", size=14, face="bold", hjust = 0.5),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))
plot(p)

# Heatmap draw
library(gplots)
# Lean.vs.Obese.topgenes <- row.names(Lean.vs.Obese)[1:length(de.Lean.vs.Obese)]
Lean.vs.Obese.topgenes <- c(up.gene,down.gene)
i <- which(row.names(v) %in% Lean.vs.Obese.topgenes)
j <- c(which(v$targets$bmi == "Lean"), which(v$targets$bmi == "Obese"))
mycol <- colorpanel(1000,"blue","white","red")
heatmap.2(lcpm[i,j], scale="row", Rowv = TRUE, Colv=FALSE, dendrogram="none",
          key = TRUE, keysize = 1.5, density.info="none",
          labRow=row.names(v)[i], labCol=bmi[j],
          col=mycol, trace="none",
          margin=c(8,6), lhei=c(2,10), main = "Lean vs Obese Heatmap")
