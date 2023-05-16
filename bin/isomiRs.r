#!/usr/bin/env Rscript

# Command line argument processing
args = commandArgs(trailingOnly=TRUE)
if (length(args) < 4) {
    stop("Usage: isomiRs.r <mirtop rawData output> <design CSV> <FDR cutoff> <log2 fold change cutoff> <comparisons CSV(optional)>", call.=FALSE)
}
rawdata <- args[1]
design_csv <- args[2]
fdr_cutoff <- as.numeric(args[3])
lfc_cutoff <- as.numeric(args[4])

# Load / install packages
if (!require("isomiRs")) {
    install.packages("BiocManager")
    BiocManager::install("isomiRs")
    library("isomiRs")
}
if (!require("DESeq2")) {
    install.packages("BiocManager")
    BiocManager::install("DESeq2")
    library("DESeq2")
}
if (!require("pheatmap")) {
    install.packages("pheatmap")
    library("pheatmap")
}
if (!require("ggplot2")) {
    install.packages("ggplot2")
    library("ggplot2")
}
if (!require("openxlsx")) {
    install.packages("openxlsx")
    library("openxlsx")
}

# Read the design csv
grouping <- read.csv(design_csv, header=TRUE, check.names=FALSE, stringsAsFactors=FALSE)
# If group label absent, use sample label
grouping$group <- ifelse(is.na(grouping$group), grouping$sample, grouping$group)
# Filter groups so that all groups have at least 2 replicates
no.replicates <- table(grouping$group)
groups <- names(no.replicates)[no.replicates > 1]
# Load the mirtop rawData
countdata <- read.csv(rawdata, header=TRUE, sep="\t", check.names=FALSE, stringsAsFactors=FALSE)
countdata <- countdata[,c(colnames(countdata)[1:6],grouping$sample)]
countdata <- as.data.frame(countdata)

# Create coldata
coldata <- data.frame(row.names=grouping$sample, group=grouping$group)
# Load isomiRs object
ids <- IsomirDataSeqFromMirtop(countdata, coldata, pct=0.05, design=~group)
# Collect raw counts data
ids <- isoCounts(ids)
rawCounts <- counts(ids)
write.xlsx(data.frame("miRNA"=row.names(rawCounts),rawCounts,check.names=FALSE), "raw_miRNA_counts.xlsx")
# Normalize counts. Note: normalized counts are log2 values, so in some cases, need to apply 2**x
#If all samples only belong to one group, set design formula to ~1
if (length(no.replicates) == 1) {
    ids <- isoNorm(ids, formula=~1)
} else if (length(groups) > 0) {
    ids <- isoNorm(ids, formula=~group)
} else {
    ids <- isoNorm(ids)
    coldata <- NA # No need to label figures with group labels if there are no replicates
}
normCounts <- counts(ids, norm=TRUE)
write.xlsx(data.frame("miRNA"=row.names(normCounts),2**normCounts,check.names=FALSE), "normalized_miRNA_counts.xlsx")

# Make study-level figures, only happens if there are more than 2 samples
if (length(grouping$sample) > 2) {
    # Make a MDS plot of samples
    jpeg("isomiRs_sample_MDS_plot.jpg", width=8, height=8, unit="in", res=300)
    p <- limma::plotMDS(normCounts)
    dev.off()
    mdsdata <- data.frame(p$cmdscale.out)
    names(mdsdata) <- c("Dimension 1", "Dimension 2")
    mdsdata$group <- ifelse(is.na(coldata), NA, coldata$group)
    write.table(mdsdata, "isomiRs_sample_MDS_plot.tsv", sep="\t", quote=FALSE)
    # Calculate correlation matrix
    correlation_matrix <- data.frame(cor(normCounts, method="pearson"))
    jpeg("isomiRs_sample_similarity_matrix.jpg", width=8, height=8, unit="in", res=300)
    pheatmap(correlation_matrix, annotation_col=coldata)
    dev.off()
    # Also output the matrix to plot in MultiQC
    sample_order <- hclust(as.dist(1-correlation_matrix))$order
    correlation_matrix <- correlation_matrix[sample_order, sample_order]
    write.table(correlation_matrix, "isomiRs_sample_similarity_matrix.tsv", sep="\t", quote=FALSE)
    # Make a heatmap for top 100 genes with most variance
    topVarGenes <- head(order(rowVars(normCounts), decreasing=TRUE), 100)
    topGeneMat <- normCounts[topVarGenes, ]
    topGeneMat <- topGeneMat - rowMeans(topGeneMat)
    jpeg("isomiRs_top_miRNA_gene_heatmap.jpg", width=8, height=12, unit="in", res=300)
    pheatmap(topGeneMat, annotation_col=coldata, fontsize_row=8)
    dev.off()
    # Also output the heatmap matrix to plot in MultiQC
    rowclust <- hclust(dist(topGeneMat))
    colclust <- hclust(dist(t(topGeneMat)))
    topGeneMat <- topGeneMat[rowclust$order, colclust$order]
    write.table(topGeneMat, "isomiRs_top_miRNA_gene_heatmap.tsv", sep="\t", quote=FALSE)
}

# Carry out DGE analysis if there are at least 2 groups with replicates
if (length(groups)>1) {
    # Remove groups with less than 2 replicates and calculate geomean of transformed counts per group
    ids <- ids[,colData(ids)$group %in% groups]
    dds <- isoDE(ids, formula=~group)
    transMeanPerLvl <- sapply(unique(ids$group), function(lvl) rowMeans(normCounts[,ids$group==lvl]))
    transMeanPerLvl <- as.data.frame(2**transMeanPerLvl)
    # Read the group comparison file
    if (length(args) > 4) {
        comparisons_file <- args[5]
        comparisons <- read.csv(comparisons_file, header=TRUE, check.names=FALSE, stringsAsFactors=FALSE)
        comparisons <- t(comparisons)
        rownames(comparisons) <- c()
    } else {
        comparisons <- combn(groups, 2)
    }
    for (i in 1:ncol(comparisons)) {
        cond1 <- comparisons[1,i]
        cond2 <- comparisons[2,i]
        res <- results(dds, contrast=c("group",cond1,cond2), alpha=fdr_cutoff, lfcThreshold=lfc_cutoff)
        # Add mean normalized counts for each group
        res[,cond1] <- transMeanPerLvl[,cond1]
        res[,cond2] <- transMeanPerLvl[,cond2]
        # Shrink the fold change. Cannot state lfcThreshold using "ashr", need to copy padj from res to reslfc
        reslfc <- lfcShrink(dds, contrast=c("group",cond1,cond2), type="ashr")
        res$lfcShrink <- reslfc$log2FoldChange
        reslfc$padj <- res$padj
        # Sort the results by FDR and pvalue
        res <- res[order(res$padj, res$pvalue), ]
        # Save results
        write.table(res, paste0(cond1,"_vs_",cond2,"_isomiRs_results.tsv"), quote=FALSE, sep="\t")
        write.xlsx(data.frame(res), paste0(cond1,"_vs_",cond2,"_isomiRs_results.xlsx"), row.names=TRUE)
        # Save MA plot
        jpeg(paste0(cond1,"_vs_",cond2,"_isomiRs_MA_plot.jpg"), width=8, height=8, unit="in", res=300)
        plotMA(reslfc, alpha=fdr_cutoff)
        dev.off()
        # Plot scatter plot
        res <- as.data.frame(res[complete.cases(res),])
        plot_data <- res[c(cond1, cond2)]
        plot_data$DEG <- res$padj<fdr_cutoff
        jpeg(paste0(cond1,"_vs_",cond2,"_isomiRs_scatterplot.jpg"), width=8, height=8, unit="in", res=300)
        print(ggplot(plot_data, aes_string(cond1, cond2, colour="DEG"))+geom_point()+scale_x_log10(limits=c(1,NA))+scale_y_log10(limits=c(1,NA)))
        dev.off()
    }
}
