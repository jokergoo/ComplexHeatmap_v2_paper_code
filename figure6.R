library(GenomicRanges)
library(data.table)
library(EnrichedHeatmap)
library(circlize)

states_bed = fread("E003_15_coreMarks_mnemonics.bed")
states = GRanges(seqnames = states_bed[[1]], 
    ranges = IRanges(states_bed[[2]] + 1, states_bed[[3]]), 
    states = states_bed[[4]])
unique(states_bed[[4]])

map = c(
    "1_TssA"      = "TssActive",
    "2_TssAFlnk"  = "TssActive",
    "3_TxFlnk"    = "Transcript",
    "4_Tx"        = "Transcript",
    "5_TxWk"      = "Transcript",
    "6_EnhG"      = "Enhancer",
    "7_Enh"       = "Enhancer",
    "8_ZNF/Rpts"  = "Heterochromatin",
    "9_Het"       = "Heterochromatin",
    "10_TssBiv"   = "TssBivalent",
    "11_BivFlnk"  = "TssBivalent",
    "12_EnhBiv"   = "Enhancer",
    "13_ReprPC"   = "Repressive",
    "14_ReprPCWk" = "Repressive",
    "15_Quies"    = "Quiescent"
)
states$states_simplified = map[states$states]

states_col = c(
    "TssActive"       = "Red",
    "TssBivalent"     = "Orange",
    "Transcript"      = "Green",
    "Enhancer"        = "Yellow",
    "Heterochromatin" = "PaleTurquoise",
    "Repressive"      = "Grey",
    "Quiescent"       = "black"
)
states_name = names(states_col)
n_states = length(states_col)

states$states_simplified = factor(states$states_simplified, levels = states_name)


library(GenomicFeatures)
txdb = loadDb("gen10_long_protein_coding_gene_adjusted.sqlite")

g = genes(txdb)
tss = promoters(g, upstream = 0, downstream = 1)

tss_chr1 = tss[seqnames(tss) == "chr1"]

# column "states_simplified" is in character mode
mat_states = normalizeToMatrix(states, tss_chr1, value_column = "states_simplified")

expr = read.table("57epigenomes.RPKM.pc.gz", row.names = 1, header = TRUE)
expr = as.matrix(expr)
obj = readRDS("chr1_roadmap_merged_bsseq.rds")
meth = granges(obj)
meth_mat = getMeth(obj, type = "smooth")
mcols(meth) = meth_mat

names(tss_chr1) = gsub("\\.\\d+$", "", names(tss_chr1))
cn = intersect(names(tss_chr1), rownames(expr))
tss_chr1 = tss_chr1[cn]
expr = expr[cn, ]

mat_states = normalizeToMatrix(states, tss_chr1, value_column = "states_simplified")
mat_meth = normalizeToMatrix(meth, tss_chr1, value_column = "E003", mean_mode = "absolute",
    smooth = TRUE)

fi = attr(mat_meth, "failed_rows")
ind = setdiff(seq_len(length(tss_chr1)), fi)
tss_chr1 = tss_chr1[ind]
mat_states = mat_states[ind, ]
mat_meth = mat_meth[ind, ]
e = log2(expr[names(tss_chr1), "E003"] + 1)

v = rowMeans(mat_states[, 95:110])
split = ifelse(v < 1.7, "active", ifelse(v > 2.6, "inactive", "bivalent"))
split = factor(split, levels = c("active", "bivalent", "inactive"))

pdf("figure6.pdf", width = 8, height = 10)
meth_col_fun = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))
ht_list = EnrichedHeatmap(mat_states, name = "States", col = states_col, cluster_rows = TRUE, 
    row_split = split, use_raster = TRUE, raster_quality = 4,
    top_annotation = HeatmapAnnotation(enrich = anno_enriched(gp = gpar(lty = 3:1)))) +
EnrichedHeatmap(mat_meth, name = "Methylation", col = meth_col_fun, use_raster = TRUE,raster_quality = 4,
    top_annotation = HeatmapAnnotation(enrich = anno_enriched(gp = gpar(lty = 3:1)))) +
Heatmap(e, name = "Expression", use_raster = TRUE,raster_quality = 4,
    show_row_names = FALSE, width = unit(1, "cm"), show_column_names = FALSE,
    top_annotation = HeatmapAnnotation(summary = anno_summary(gp = gpar(lty = 3:1),
        outline = FALSE, axis_param = list(side = "right"))))

lgd = Legend(title = "Category", legend_gp = gpar(lty = 1:3), grid_width = unit(1, "cm"),
    labels = levels(split), type = "lines")
draw(ht_list, ht_gap = unit(8, "mm"), column_title = "Visualize chromatin states, methylation and expression",
    heatmap_legend_list = list(lgd))
dev.off()

