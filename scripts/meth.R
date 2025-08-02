library(TCGAbiolinks)
library(here)
library(knitr)
library(limma)
library(minfi)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylation450kmanifest)
library(RColorBrewer)
library(missMethyl)
library(minfiData)
library(Gviz)
library(DMRcate)
library(stringr)
library(survival)
library(survminer)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(sesame)
library(sesameData)
library(SummarizedExperiment)
library(ggplot2)
library(GenomicRanges)
library(pheatmap)
library(manipulate)
library(ggpubr)
library(gghalves)
library(plotly)
library(dplyr)
library(ggplot2)
library(reshape2)
# retrieve data 
here()
query_met <- GDCquery(
  project = "TARGET-WT",
  data.category = "DNA Methylation",
  data.type = "Methylation Beta Value",
  platform = "Illumina Human Methylation 450"
)
#  Download the data
GDCdownload(query_met, method = "api", files.per.chunk = 5)
#  Prepare the data
met_data <- GDCprepare(query_met)
# save
beta_mat <- assay(met_data)
coldata <- colData(met_data)
metadata <- as.data.frame(coldata)
saveRDS(met_data, file = "data/TCGA_met_data.rds")

# assign groups based on our hypothesis 

metadata <- metadata %>%
  mutate(
    ten_year_survival = case_when(
      is.na(vital_status) | (vital_status == "Dead" & is.na(days_to_death)) ~ NA_character_,
      vital_status == "Dead" & days_to_death <= 10*365 ~ "dead",
      vital_status == "Dead" & days_to_death > 10*365 ~ "alive",
      vital_status == "Alive" ~ "alive",
      TRUE ~ NA_character_
    ),
    ten_year_survival = factor(ten_year_survival, levels = c("alive", "dead"))
  )
metadata$ten_year_survival
ggplot(metadata, aes(x = ten_year_survival)) +
  geom_bar(position = "dodge") +
  labs(
    title = "ten_year_survival",
    x = "ten_year_survival",
    y = "Number of Patients"
  ) +
  theme_minimal()
sum(is.na(metadata$ten_year_survival))
metadata <- metadata[!is.na(metadata$ten_year_survival),]
beta_mat <- beta_mat[, colnames(beta_mat) %in% metadata$barcode]

# quality control
densityPlot(beta_mat, sampGroups = metadata$ten_year_survival) #density plot to examine beta values
plotMDS(beta_mat, top = 1000, labels = metadata$ten_year_survival, 
  col = as.numeric(as.factor(metadata$ten_year_survival)), gene.selection = "common", main = "MDS plot of TARGET-WT samples")
legend("topright", legend = levels(as.factor(metadata$ten_year_survival)),
 col = 1:length(unique(metadata$ten_year_survival)), pch = 16)

## Filter probes
# Load 450k annotationann450k = getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann450k = getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
# Remove probes with any NA
table(rowSums(is.na(beta_mat))==0) # examine number of probes with at one or more NA
beta_mat = beta_mat[rowSums(is.na(beta_mat)) == 0, ]
# Remove cross-reactive probes (Chen et al. 2013)
# Download and remove cross-reactive probes
url <- "https://raw.githubusercontent.com/hamidghaedi/Methylation_Analysis/master/cross_reactive_probe.chen2013.csv"
download.file(url, destfile = "cross_reactive_probe.chen2013.csv", mode = "wb")
cross_reactive = read.csv("cross_reactive_probe.chen2013.csv")
cross_reactive_ids = cross_reactive$TargetID[-1]  # remove header
beta_mat = beta_mat[!rownames(beta_mat) %in% cross_reactive_ids, ]
## QC plots with filtered beta values matrix
densityPlot(beta_mat, sampGroups = metadata$ten_year_survival) #density plot to examine beta values
plotMDS(beta_mat, top = 1000, labels = metadata$ten_year_survival, 
  col = as.numeric(as.factor(metadata$ten_year_survival)), gene.selection = "common", main = "MDS plot of TARGET-WT samples")
legend("topright", legend = levels(as.factor(metadata$ten_year_survival)),
 col = 1:length(unique(metadata$ten_year_survival)), pch = 16)
# Transforms beta values to M values using the logit transformation: M = log2(β / (1 - β)) with an offset to prevent infinite values when β = 0 or 1.
beta2m = function(beta_mat, offset = 1e-6) {
  beta_mat = pmin(pmax(beta_mat, offset), 1 - offset)
  log2(beta_mat / (1 - beta_mat))
}

mval_mat = beta2m(beta_mat)

# Or 
# mval_mat = BetaValueToMValue(beta_mat)
# Save Cleaned Matrices
saveRDS(beta_mat, "TCGA-WT_filtered.rds")
saveRDS(mval_mat, "TCGA-WT_mvalue_filtered.rds")
#

# Differential methylation analysis
# remove the patients dropped 
colnames(mval_mat)
metadata$barcode
mval_mat <- mval_mat[, colnames(mval_mat) %in% metadata$barcode]
# 
design = model.matrix(~ metadata$ten_year_survival)
fit = lmFit(mval_mat, design)
fit2 = eBayes(fit)

# extract results with annotation
ann450k_sub <- ann450k[match(rownames(mval_mat), ann450k$Name), 
                       c(1:4, 12:19, 24:ncol(ann450k))]
DMPs = topTable(fit2, coef=2, number=Inf, genelist = ann450k_sub)
head(DMPs)

# plot top 10 differentially methylated probes
par(mfrow=c(2,5))
sapply(rownames(DMPs)[1:10], function(cpg){
plotCpg(beta_mat, cpg=cpg, pheno= metadata$ten_year_survival, ylab="Beta values")
})

#
## Differentially methylated regions
my_annotation = cpg.annotate(object = mval_mat, datatype = "array", what = "M", analysis.type = "differential", design = design, contrasts = FALSE, coef = 2, arraytype = "450K")
my_annotation
DMRs = dmrcate(my_annotation, lambda = 1000, C = 2)
result.ranges = extractRanges(DMRs)
result.ranges[1]$overlapping.genes
#

# Convert result.ranges to a data.frame
dmr_df <- as.data.frame(result.ranges)
head(dmr_df)
# Rank by FDR (Stouffer p-value)
dmr_df_ranked <- dmr_df[order(dmr_df$Stouffer, decreasing = FALSE), ]
# View top 5 ranked DMRs
head(dmr_df_ranked, 5)
# DMRs by chromosome
table(seqnames(result.ranges))


# vis
### custom plot of DMRs using Gviz
gen = "hg19" #genome version to be used
dmrIndex = 1 # will use ch with genes GATA3
groups <- metadata$ten_year_survival 
groups <- as.factor(groups)
# or any other grouping variable
length(groups)  # should now be 130
levels(groups)


# sample group colors
pal = brewer.pal(n = length(levels(groups)), name = "Set1")
names(pal) = levels(groups)
cols = pal[groups]


# Extract region of interest
chrom = as.character(seqnames(result.ranges[dmrIndex]))
start = as.numeric(start(result.ranges[dmrIndex]))
end   = as.numeric(end(result.ranges[dmrIndex]))


# Add some padding (25% extra space) to view context
minbase = start - 0.25 * (end - start)
maxbase = end   + 0.25 * (end - start)

## load annotation tracks
# CpG islands file (from UCSC)
islands = read.table("cpgIslandExt.txt", header=FALSE, stringsAsFactors=FALSE)
islandData = GRanges(seqnames = islands[,2],
              ranges = IRanges(start=islands[,3], end=islands[,4]),
              strand = "*")


# filter CpG islands to region of interest
roi <- GRanges(seqnames = chrom, ranges = IRanges(start = minbase, end = maxbase))
islandData_sub <- subsetByOverlaps(islandData, roi)

# DNase hypersensitive sites file (from UCSC)
dnase_df = read.table("wgEncodeRegDnaseClusteredV3.bed", header=FALSE, stringsAsFactors=FALSE)

# Name the important columns
colnames(dnase_df)[1:3] <- c("chr", "start", "end")
colnames(dnase_df)[5] <- "score"

# Convert to GRanges
dnaseData <- GRanges(
  seqnames = dnase_df$chr,
  ranges = IRanges(start = dnase_df$start, end = dnase_df$end),
  score = as.numeric(dnase_df$score)
)

# filter DNase data to the region of interest
roi <- GRanges(seqnames = chrom, ranges = IRanges(start = minbase, end = maxbase))
dnaseData_sub <- subsetByOverlaps(dnaseData, roi)

## prepare methylation data

# make sure annotation and beta matrix are in same order
ann450kOrd = ann450k[order(ann450k$chr, ann450k$pos), ]
bValsOrd  = beta_mat[match(ann450kOrd$Name, rownames(beta_mat)), ]

# extract probes overlapping the DMR
cpgData = GRanges(seqnames = ann450kOrd$chr,ranges = IRanges(start=ann450kOrd$pos, end=ann450kOrd$pos), strand = "*", betas = bValsOrd)
cpgData = subsetByOverlaps(cpgData, result.ranges[dmrIndex])

## Create Gviz tracks
# ideogram and axis
iTrack = IdeogramTrack(genome=gen, chromosome=chrom, name="")
gTrack = GenomeAxisTrack(col="black", cex=1, name="", fontcolor="black")

# RefSeq track
rTrack = UcscTrack(genome=gen, chromosome=chrom, track ="NCBI RefSeq", table = "refGene", from=minbase, to=maxbase,  trackType="GeneRegionTrack", rstarts="exonStarts", rends="exonEnds", gene="name", symbol="name2", transcript="name", strand="strand", fill="darkblue", stacking="squish", name="RefSeq", showId=TRUE, geneSymbol=TRUE)

# CpG islands track
islandTrack = AnnotationTrack(range=islandData_sub, genome=gen, name="CpG Is.", chromosome=chrom, fill="darkgreen")

# DNase hypersensitive sites track
dnaseTrack = AnnotationTrack(range=dnaseData_sub, genome=gen, name="DNaseI", chromosome=chrom, fill="orange")

# DMR track
dmrTrack = AnnotationTrack(start=start, end=end, genome=gen, name="DMR", chromosome=chrom, fill="red")
# Methylation data track


methTrack = DataTrack(range=cpgData, genome=gen, chromosome=chrom,
                       groups=groups, type=c("smooth"),
                       col=pal, name="Beta values", legend=TRUE,
                       background.panel="white", ylim=c(-0.05,1.05),
                       cex.title=0.8, cex.axis=0.8, cex.legend=0.8)
# type can be a, p, or smooth (better)
## combine all tracks and plot
tracks = list(iTrack, gTrack, methTrack, dmrTrack, islandTrack, dnaseTrack, rTrack)
sizes  = c(2, 2, 5, 2, 2, 2, 3)  # relative heights

plotTracks(tracks, from=minbase,
  to=maxbase, showTitle=TRUE, add53=TRUE, add35=TRUE,lty.grid=3, sizes=sizes, main="DMR1 identified on Chromosome 10")

#
dmr_beta <- beta_mat[cpgs, ]
group_means <- apply(dmr_beta, 2, function(x) mean(x, na.rm = TRUE))
boxplot(group_means ~ groups,
        ylab = "Average Beta across DMR",
        col = c("skyblue", "tomato"),
        main = paste("Mean Beta for DMR", dmrIndex))

## GO analysis
# get significant probes with adjusted p value <0.05
sigCpGs <- DMPs$Name[DMPs$adj.P.Val<0.05]
sigCpGs[1:10]

# get all probes 
all = DMPs$Name

# run enrichment
gst <- gometh(sig.cpg=sigCpGs, all.cpg=all)
gsa <- topGSA(gst)
top_terms <- gsa[1:10, ]  # top 10 enriched GO terms

ggplot(top_terms, aes(x = reorder(TERM, -P.DE), y = -log10(P.DE))) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(
    title = "Top 10 GO Terms from gometh",
    x = "GO Term",
    y = "-log10(p-value)"
  ) +
  theme_minimal()
