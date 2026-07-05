suppressWarnings(suppressMessages(library(DESeq2)))
cm <- read.csv("TICTOC_run1_filteredCounts_v3.csv", row.names=1, check.names=FALSE)
cts <- as.matrix(cm); mode(cts) <- "integer"
grp <- sub("(__r[0-9]+|[0-9]+)$","",colnames(cts)); colnames(cts) <- make.unique(colnames(cts), sep="__r")
cd <- data.frame(row.names=colnames(cts),
                 genotype=sub("(Root|Shoot).*","",grp), tissue=sub(".*?(Root|Shoot).*","\1",grp),
                 treatment=sub(".*(Flight|Ground)$","\1",grp))
dds <- DESeqDataSetFromMatrix(cts, cd, ~1); vsd <- vst(dds, blind=TRUE)
pc <- plotPCA(vsd, intgroup=c("tissue","treatment","genotype"), returnData=TRUE)
pv <- round(100*attr(pc,"percentVar"))
write.csv(pc, "deseq2/pca_coords.csv", row.names=FALSE)
pdf("manuscript/Fig3_PCA.pdf", width=7, height=5.5)
cols <- c(Root="#1b7837", Shoot="#762a83"); pch <- c(Flight=19, Ground=1)
plot(pc$PC1, pc$PC2, col=cols[pc$tissue], pch=pch[pc$treatment], cex=1.6, lwd=2,
     xlab=paste0("PC1: ",pv[1],"% variance"), ylab=paste0("PC2: ",pv[2],"% variance"),
     main="TICTOC RNA-seq PCA (VST)")
text(pc$PC1, pc$PC2, labels=pc$genotype, pos=3, cex=0.5, col="grey30")
legend("topright", c("Root","Shoot","Flight (filled)","Ground (open)"),
       col=c(cols,"black","black"), pch=c(15,15,19,1), bty="n", cex=0.9)
dev.off()
cat("PC1:",pv[1],"% PC2:",pv[2],"%\n")
