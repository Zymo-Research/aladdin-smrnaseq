#!/usr/bin/env Rscript

#Command line arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args) <3) {
	stop("Usage: tximport_deseq.r <rsem data output> <design CSV> <tx2gene>", call.=FALSE)

}
rsemquant <- args[1]
design_csv <- args[2]
fdr_cutoff <- as.numeric(args[3])
lfc_cutoff <- as.numeric(args[4])
tx2gene <- args[5]

#Load packages
library("dplyr")
library("tximport")
library("DESeq2")
library("pheatmap")
library("ggplot2")
library("openxlsx")

grouping <- read.csv(design_csv, header=TRUE, check.names=FALSE, stringsAsFactors = FALSE)
#If group label absent, use sample label
grouping$group <- ifelse(is.na(grouping$group), grouping$sample, grouping$group)
grouping <- grouping[,c("sample", "group")]
#Arrange the order of the samples; for deseq2 the sample order of the design file must match the sample order of the tximport object 
grouping <- grouping[order(grouping$sample),]
rownames(grouping) <- grouping$sample
#Create grouplist and remove replicates from the list
no.replicates <- table(grouping$group)
groups <- names(no.replicates)[no.replicates > 1]

#Load rsem data
#Aggregate rsem log files together into table with tximport
rsemquant <- as.data.frame(read.table(rsemquant, header=F))
#remove file path and keep filenames
rsemquant$rownames <- gsub(".*(?=\\/).", "", rsemquant$V1, perl=TRUE)
rsemquant <- rsemquant[order(rsemquant$rownames),]
rownames(rsemquant) <- rsemquant$rownames
rsemquant <- subset(rsemquant, select = -c(rownames))
rsemquant <- t(rsemquant)
#Convert dataframe to character vector
rsemquant <- rsemquant[1,]

#Use tximport to compile rsem quantification results
txi.rsem.matrix <- tximport(rsemquant, type = "rsem", txOut = TRUE, countsFromAbundance = "lengthScaledTPM")
#rsem estimated counts will be used to estimate RNA Type counts
write.table(txi.rsem.matrix$counts, "rsem_estimated_transcriptcounts.tsv", sep = "\t")
#Summarize transcripts to gene level
tx2gene <- read.table(tx2gene, sep="\t")

txi.rsem.matrix <- summarizeToGene(txi.rsem.matrix, tx2gene[,c("transcript_id", "gene_id")], ignoreTxVersion=TRUE, countsFromAbundance = "lengthScaledTPM")

#If all samples belong to one group, set design to ~1
if (length(no.replicates)==1) {
    ids <- DESeqDataSetFromTximport(txi.rsem.matrix, colData = grouping, design = ~1)
} else {
    ids <- DESeqDataSetFromTximport(txi.rsem.matrix, colData = grouping, design = ~group)
}

#Clean tx2gene to bind to rawCounts and normCounts
tx2gene <- dplyr::distinct(tx2gene, gene_id, .keep_all=T)
rownames(tx2gene) <- tx2gene$gene_id
tx2gene <- dplyr::select(tx2gene, -c(transcript_id, gene_id))

#Merge rawCounts with gene names; because rawCounts will not be used downstream, character vectors gene_id and gene_name can be added
rawCounts <- counts(ids)
rawCounts <- merge(tx2gene, rawCounts, by = 0)
colnames(rawCounts)[1] <- "gene_id"
#Raw counts saved as downloadable xlsx
write.xlsx(data.frame(rawCounts,check.names=FALSE), "rawgenes_exceptmiRNArRNA_counts.xlsx")

#Calculate normalized counts
#If no groups have more than one sample, rlog must be blind to groups or will throw an error
if (length(groups) > 0) {
    normCounts <- assay(rlog(ids, blind = F))
} else {
    normCounts <- assay(rlog(ids, blind = T))
}

normCounts_xlsx <- data.frame(2**normCounts, check.names=FALSE)
normCounts_xlsx <- merge(tx2gene, normCounts_xlsx, by = 0)
colnames(normCounts_xlsx)[1] <- "gene_id"
write.xlsx(normCounts_xlsx, "normalizedgenes_exceptmiRNArRNA_counts.xlsx")
rm(normCounts_xlsx)

if (length(grouping$sample) > 2) {
    # Make a MDS plot of samples
    jpeg("deseq2_sample_MDS_plot.jpg", width=8, height=8, unit="in", res=300)
    p <- limma::plotMDS(normCounts)
    dev.off()

    mdsdata <- data.frame(p[c("x", "y")])
    colnames(mdsdata) <- c("Dimension 1", "Dimension 2")
    mdsdata$group <- ifelse(is.na(grouping$group), NA, grouping$group)
    write.table(mdsdata, "deseq2_sample_mds_plot.tsv", sep="\t", quote=FALSE)

    #Calculate correlation matrix
    correlation_matrix <- data.frame(cor(normCounts, method="pearson"))
    jpeg("deseq2_sample_similarity_matrix.jpg", width=8, height=8, unit="in", res=300)
    pheatmap(correlation_matrix, annotation_col=grouping)
    dev.off()
    sample_order <- hclust(as.dist(1-correlation_matrix))$order
    correlation_matrix <- correlation_matrix[sample_order, sample_order]
    write.table(correlation_matrix, "deseq2_sample_similarity_matrix.tsv", sep="\t", quote=FALSE)
    #Make a heatmap for top 100 genes with the most variance
    topVarGenes <- head(order(rowVars(normCounts), decreasing=TRUE), 100)
    topGeneMat <- normCounts[topVarGenes, ]
    topGeneMat <- topGeneMat - rowMeans(topGeneMat)
    jpeg("deseq2_top_nonmiRNA_gene_heatmap.jpg", width=8, height=12, unit="in", res=300)
    pheatmap(topGeneMat, annotation_col=grouping, fontsize_row=8)
    dev.off()
    #Also output the heatmap matrix to plot in MultiQC
    rowclust <- hclust(dist(topGeneMat))
    colclust <- hclust(dist(t(topGeneMat)))
    topGeneMat <- topGeneMat[, colclust$order]
    order <- rownames(topGeneMat[rowclust$order,])
    topGeneMat <- merge(tx2gene, topGeneMat, by=0)
    rownames(topGeneMat) <- topGeneMat$Row.names
    topGeneMat <- select(topGeneMat, -c(Row.names))
    topGeneMat <- topGeneMat[order,] 
    write.table(topGeneMat, "deseq2_top_nonmiRNA_gene_heatmap.tsv", sep="\t", quote=FALSE)
}

#Generate results table for use in deseq2 module
if (length(groups) > 1){
    #Remove groups with less than 2 replicates and calculate geomean of transformed counts per group
    ids <- ids[,ids$group %in% groups]
    ids$group <- droplevels(ids$group)
    ids$group <- factor(ids$group, levels = groups)
    dds <- DESeq(ids)
    transMeanPerLvl <- sapply(unique(ids$group), function(lvl) rowMeans(normCounts[,ids$group==lvl]))
    transMeanPerLvl <- as.data.frame(2**transMeanPerLvl)
    colnames(transMeanPerLvl) <- levels(unique(ids$group))
    #Read the group comparison file
    if (length(args) > 5){
        comparisons_file <- args[6]
        comparisons <- read.csv(comparisons_file, header=TRUE, check.names=FALSE, stringsAsFactors=FALSE)
        comparisons <- t(comparisons)
        rownames(comparisons) <- c()
    } else {
        comparisons <- combn(groups, 2)
    }
    for(i in 1:ncol(comparisons)){
        cond1 <- comparisons[1, i]
        cond2 <- comparisons[2, i]
        res <- results(dds, contrast=c("group",cond1, cond2), alpha=fdr_cutoff, lfcThreshold=lfc_cutoff)
        #Add mean normalized counts for each group
        res[,cond1] <- transMeanPerLvl[,cond1]
        res[,cond2] <- transMeanPerLvl[,cond2]
        #Shrink the fold change. Cannot state lfcThreshold using "ashr", need to copy padj from res to reslfc
        reslfc <- lfcShrink(dds, contrast=c("group",cond1,cond2), type="ashr")
        res$lfcShrink <- reslfc$log2FoldChange
        #Line below removed; was causing problems
        reslfc$padj <- res$padj
	#Prepare res for printing into deseq2 results comparison file
	res <- merge(tx2gene, data.frame(res), by = 0)
        #Sort the results by FDR and pvalue
        res <- res[order(res$padj, res$pvalue),]
        #Save results
	rownames(res) <- res$Row.names
	res <- select(res, -c("Row.names"))
        write.table(res, paste0(cond1,"_vs_",cond2,"_deseq2_results.tsv"), quote=FALSE, sep="\t")
	write.xlsx(res, paste0(cond1,"_vs_",cond2,"_deseq2_results.xlsx"), row.names=TRUE)
        #Save MA plot
        jpeg(paste0(cond1,"_vs_",cond2,"_deseq2_MA_plot.jpg"), width=8, height=8, unit="in", res=300)
        plotMA(reslfc, alpha=fdr_cutoff)
        dev.off()
        #plot scatter plot
        res <- as.data.frame(res[complete.cases(res),])
        plot_data <- res[c(cond1, cond2)]
        plot_data$DEG <- res$padj<fdr_cutoff
        jpeg(paste0(cond1,"_vs_",cond2,"_deseq2_scatterplot.jpg"), width=8, height=8, unit="in", res=300)
        print(ggplot(plot_data, aes_string(cond1, cond2, colour="DEG"))+geom_point()+scale_x_log10(limits=c(1,NA))+scale_y_log10(limits=c(1,NA)))
        dev.off()
    }
}
