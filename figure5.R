
library(ComplexHeatmap)

set.seed(123)

m = cbind(matrix(rnorm(10*100), ncol = 10),
          matrix(runif(10*100, min = -2, max = 2) + 0.5, ncol = 10))
colnames(m) = paste0("C", 1:ncol(m))

ha1 = HeatmapAnnotation(distribution = c(rep("rnorm", 10), rep("runif", 10)), col = list(distribution = c("rnorm" = 2, "runif" =3)))
p1 = grid.grabExpr(draw(densityHeatmap(m, ylab = "Value", top_annotation = ha1, column_title = "A) Density heatmap"), merge_legends = TRUE))


p2 = grid.grabExpr(draw(frequencyHeatmap(m, ylab = "Value", top_annotation = ha1, use_3d = TRUE, column_title = "B) 3D frequency heatmap")))


mat = read.table(system.file("extdata", package = "ComplexHeatmap", 
    "tcga_lung_adenocarcinoma_provisional_ras_raf_mek_jnk_signalling.txt"), 
    header = TRUE, stringsAsFactors = FALSE, sep = "\t")
mat[is.na(mat)] = ""
rownames(mat) = mat[, 1]
mat = mat[, -1]
mat=  mat[, -ncol(mat)]
mat = t(as.matrix(mat))
mat = mat[, sample(ncol(mat), 80)]

col = c("HOMDEL" = "blue", "AMP" = "red", "MUT" = "#008000")
alter_fun = list(
    background = function(x, y, w, h) {
        grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
            gp = gpar(fill = "#CCCCCC", col = NA))
    },
    # big blue
    HOMDEL = function(x, y, w, h) {
        grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
            gp = gpar(fill = col["HOMDEL"], col = NA))
    },
    # big red
    AMP = function(x, y, w, h) {
        grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
            gp = gpar(fill = col["AMP"], col = NA))
    },
    # small green
    MUT = function(x, y, w, h) {
        grid.rect(x, y, w-unit(2, "pt"), h*0.33, 
            gp = gpar(fill = col["MUT"], col = NA))
    }
)

column_title = "C) OncoPrint"
heatmap_legend_param = list(title = "Alternations", at = c("HOMDEL", "AMP", "MUT"), 
        labels = c("Deep deletion", "Amplification", "Mutation"), title_position = "leftcenter", nrow = 1)
sample_order = scan(paste0(system.file("extdata", package = "ComplexHeatmap"), 
    "/sample_order.txt"), what = "character")
p3 = grid.grabExpr(draw(oncoPrint(mat,
    alter_fun = alter_fun, col = col, 
    row_order = 1:nrow(mat), column_order = sample_order,
    remove_empty_columns = TRUE, remove_empty_rows = TRUE,
    column_title = column_title, heatmap_legend_param = heatmap_legend_param), heatmap_legend_side = "bottom"))




file_list = c(
    "ESC" = "~/project/development/ComplexHeatmap-reference/data/E016-H3K4me3.narrowPeak.gz",
    "ES-deriv1" = "~/project/development/ComplexHeatmap-reference/data/E004-H3K4me3.narrowPeak.gz",
    "ES-deriv2" = "~/project/development/ComplexHeatmap-reference/data/E006-H3K4me3.narrowPeak.gz",
    "Brain" = "~/project/development/ComplexHeatmap-reference/data/E071-H3K4me3.narrowPeak.gz",
    "Muscle" = "~/project/development/ComplexHeatmap-reference/data/E100-H3K4me3.narrowPeak.gz",
    "Heart" = "~/project/development/ComplexHeatmap-reference/data/E104-H3K4me3.narrowPeak.gz"
)
library(GenomicRanges)
peak_list = lapply(file_list, function(f) {
    df = read.table(f)
    GRanges(seqnames = df[, 1], ranges = IRanges(df[, 2], df [, 3]))
})

m = make_comb_mat(peak_list)
m = m[comb_size(m) > 500000]

p4 = grid.grabExpr(draw(UpSet(m, column_title = "D) UpSet plot")))



library(circlize)
library(GenomicRanges)
chr_df = read.chromInfo()$df
chr_df = chr_df[chr_df$chr %in% paste0("chr", 1:22), ]
chr_gr = GRanges(seqnames = chr_df[, 1], ranges = IRanges(chr_df[, 2] + 1, chr_df[, 3]))

library(EnrichedHeatmap)
chr_window = makeWindows(chr_gr, w = 1e6)

average_in_window = function(window, gr, v, method = "weighted", empty_v = NA) {

    if(missing(v)) v = rep(1, length(gr))
    if(is.null(v)) v = rep(1, length(gr))
    if(is.atomic(v) && is.vector(v)) v = cbind(v)

    v = as.matrix(v)
    if(is.character(v) && ncol(v) > 1) {
        stop("`v` can only be a character vector.")
    }

    if(length(empty_v) == 1) {
        empty_v = rep(empty_v, ncol(v))
    }

    u = matrix(rep(empty_v, each = length(window)), nrow = length(window), ncol = ncol(v))

    mtch = as.matrix(findOverlaps(window, gr))
    intersect = pintersect(window[mtch[,1]], gr[mtch[,2]])
    w = width(intersect)
    v = v[mtch[,2], , drop = FALSE]
    n = nrow(v)

    ind_list = split(seq_len(n), mtch[, 1])
    window_index = as.numeric(names(ind_list))
    window_w = width(window)

    if(is.character(v)) {
        for(i in seq_along(ind_list)) {
            ind = ind_list[[i]]
            if(is.function(method)) {
                u[window_index[i], ] = method(v[ind], w[ind], window_w[i])
            } else {
                tb = tapply(w[ind], v[ind], sum)
                u[window_index[i], ] = names(tb[which.max(tb)])
            }
        }
    } else {
        if(method == "w0") {
            gr2 = reduce(gr, min.gapwidth = 0)
            mtch2 = as.matrix(findOverlaps(window, gr2))
            intersect2 = pintersect(window[mtch2[, 1]], gr2[mtch2[, 2]])

            width_intersect = tapply(width(intersect2), mtch2[, 1], sum)
            ind = unique(mtch2[, 1])
            width_setdiff = width(window[ind]) - width_intersect

            w2 = width(window[ind])

            for(i in seq_along(ind_list)) {
                ind = ind_list[[i]]
                x = colSums(v[ind, , drop = FALSE]*w[ind])/sum(w[ind])
                u[window_index[i], ] = (x*width_intersect[i] + empty_v*width_setdiff[i])/w2[i]
            }

        } else if(method == "absolute") {
            for(i in seq_along(ind_list)) {
                u[window_index[i], ] = colMeans(v[ind_list[[i]], , drop = FALSE])
            }
            
        } else if(method == "weighted") {
            for(i in seq_along(ind_list)) {
                ind = ind_list[[i]]
                u[window_index[i], ] = colSums(v[ind, , drop = FALSE]*w[ind])/sum(w[ind])
            }
        } else {
            if(is.function(method)) {
                for(i in seq_along(ind_list)) {
                    ind = ind_list[[i]]
                    u[window_index[i], ] = method(v[ind], w[ind], window_w[i])
                }
            } else {
                stop("wrong method.")
            }
        }
    }

    return(u)
}

bed1 = generateRandomBed(nr = 1000, nc = 10) # generateRandomBed() is from circlize package
# convert to a GRanes object
gr1 = GRanges(seqnames = bed1[, 1], ranges = IRanges(bed1[, 2], bed1[, 3]))

num_mat = average_in_window(chr_window, gr1, bed1[, -(1:3)])


bed_list = lapply(1:10, function(i) {
    generateRandomBed(nr = 1000, nc = 1, 
        fun = function(n) sample(c("gain", "loss"), n, replace = TRUE))
})
char_mat = NULL
for(i in 1:10) {
    bed = bed_list[[i]]
    bed = bed[sample(nrow(bed), 20), , drop = FALSE]
    gr_cnv = GRanges(seqnames = bed[, 1], ranges = IRanges(bed[, 2], bed[, 3]))

    char_mat = cbind(char_mat, average_in_window(chr_window, gr_cnv, bed[, 4]))
}

bed2 = generateRandomBed(nr = 100, nc = 2)
gr2 = GRanges(seqnames = bed2[, 1], ranges = IRanges(bed2[, 2], bed2[, 3]))

v = average_in_window(chr_window, gr2, bed2[, 4:5])

bed3 = generateRandomBed(nr = 40, nc = 0)
gr3 = GRanges(seqnames = bed3[, 1], ranges = IRanges(bed3[, 2], bed3[, 2]))
gr3$gene = paste0("gene_", 1:length(gr3))

mtch = as.matrix(findOverlaps(chr_window, gr3))
at = mtch[, 1]
labels = mcols(gr3)[mtch[, 2], 1]

chr = as.vector(seqnames(chr_window))
chr_level = paste0("chr", 1:22)
chr = factor(chr, levels = chr_level)

subgroup = rep(c("A", "B"), each = 5)

ht_opt$TITLE_PADDING = unit(c(4, 4), "points")
ht_list = Heatmap(t(num_mat), name = "mat", col = colorRamp2(c(-1, 0, 1), c("green", "white", "red")),
    column_split = chr, cluster_columns = FALSE, show_row_dend = FALSE,
    row_split = subgroup, cluster_row_slices = FALSE, 
    left_annotation = rowAnnotation(subgroup = subgroup, show_annotation_name = FALSE,
        annotation_legend_param = list(
            subgroup = list(direction = "horizontal", title_position = "lefttop", nrow = 1))),
    column_title_gp = gpar(fontsize = 10), border = TRUE,
    column_gap = unit(0, "points"),
    column_title = ifelse(1:22 %% 2 == 0, paste0("\n", chr_level), paste0(chr_level, "\n")),
    height = unit(1, "null"),
    heatmap_legend_param = list(direction = "horizontal", title_position = "lefttop")) %v%
Heatmap(t(char_mat), name = "CNV", col = c("gain" = "red", "loss" = "blue"),
    border = TRUE, na_col = "#EFEFEF",
    row_title = "CNV",
    height = unit(0.6, "null"),
    heatmap_legend_param = list(direction = "horizontal", title_position = "lefttop", nrow = 1)) %v%
HeatmapAnnotation(label = anno_mark(at = at, labels = labels, side = "bottom")) %v%
HeatmapAnnotation(points = anno_points(v, gp = gpar(col = 4:5), pch = c(1, 16)),
    annotation_name_side = "left", height = unit(2, "cm")) %v%
HeatmapAnnotation(bars = anno_barplot(v[, 1], gp = gpar(col = ifelse(v[ ,1] > 0, 2, 3))), 
    annotation_name_side = "left", height = unit(2, "cm"))
p5 = grid.grabExpr(draw(ht_list, heatmap_legend_side = "bottom", merge_legend = TRUE,
    column_title = "E) Genome-level multi-track heatmap"), width = 14)



library(cowplot)

pdf("figure5.pdf", width = 12, height = 15)
p = plot_grid( plot_grid(p1, p2, nrow = 1),
            plot_grid(p3, p4, nrow = 1),
            p5, nrow = 3, rel_heights = c(1, 1.2, 1.8))
print(p)
dev.off()

