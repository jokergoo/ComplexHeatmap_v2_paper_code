library(ComplexHeatmap)
library(circlize)
library(cowplot)
library(GetoptLong)

### figure A

set.seed(123)
nr1 = 4; nr2 = 8; nr3 = 6; nr = nr1 + nr2 + nr3
nc1 = 6; nc2 = 8; nc3 = 10; nc = nc1 + nc2 + nc3
mat = cbind(rbind(matrix(rnorm(nr1*nc1, mean = 1,   sd = 0.5), nr = nr1),
          matrix(rnorm(nr2*nc1, mean = 0,   sd = 0.5), nr = nr2),
          matrix(rnorm(nr3*nc1, mean = 0,   sd = 0.5), nr = nr3)),
    rbind(matrix(rnorm(nr1*nc2, mean = 0,   sd = 0.5), nr = nr1),
          matrix(rnorm(nr2*nc2, mean = 1,   sd = 0.5), nr = nr2),
          matrix(rnorm(nr3*nc2, mean = 0,   sd = 0.5), nr = nr3)),
    rbind(matrix(rnorm(nr1*nc3, mean = 0.5, sd = 0.5), nr = nr1),
          matrix(rnorm(nr2*nc3, mean = 0.5, sd = 0.5), nr = nr2),
          matrix(rnorm(nr3*nc3, mean = 1,   sd = 0.5), nr = nr3))
   )
mat = mat[sample(nr, nr), sample(nc, nc)] # random shuffle rows and columns
rownames(mat) = paste0("row", seq_len(nr))
colnames(mat) = paste0("column", seq_len(nc))

library(dendextend)
column_dend = as.dendrogram(hclust(dist(t(mat))))
column_dend = color_branches(column_dend, k = 3) # `color_branches()` returns a dendrogram object
column_dend = dendrapply(column_dend, function(d) {
    if(runif(1) > 0.5) attr(d, "nodePar") = list(cex = 0.5, pch = sample(20, 1), col = rand_color(1))
    return(d)
})
ht = Heatmap(mat, name = "mat", row_dend_width = unit(2, "cm"), cluster_columns = column_dend,
    column_title = "A) A heatmap with various annotations",
    show_column_names = FALSE,
	column_names_rot = 45, row_split = rep(c("A", "B"), 9), row_km = 2, column_split = 3,
	top_annotation = HeatmapAnnotation(foo1 = 1:24, bar1 = anno_points(runif(24))),
    right_annotation = rowAnnotation(foo2 = 18:1, bar2 = anno_barplot(cbind(runif(18), runif(18)), gp = gpar(fill = 2:3), width = unit(2, "cm"))))
p1 = grid.grabExpr({
	draw(ht, annotation_legend_list = list(Legend(title = "bar2", labels = c("group1", "group2"), legend_gp = gpar(fill = 2:3))),
		merge_legend = TRUE, padding = unit(c(5, 5, 5, 5), "mm"))
	for(slice in 1:3) {
		decorate_annotation("bar1", slice = slice, {
			grid.lines(c(0, 1), c(0.5, 0.5), gp = gpar(lty = 2, col = "#AAAAAA"))
		})
	}
})

### figure B

small_mat = mat[1:9, 1:9]
col_fun = colorRamp2(c(-2, 0, 2), c("green", "white", "red"))
ht = Heatmap(small_mat, name = "mat", col = col_fun,
    column_title = "B) Customize heatmap body",
    row_km = 2, column_km = 2,
    layer_fun = function(j, i, x, y, w, h, fill, slice_r, slice_c) {
        # restore_matrix() is explained after this chunk of code
        ind_mat = restore_matrix(j, i, x, y)
        for(ir in seq_len(nrow(ind_mat))) {
            # start from the second column
            for(ic in seq_len(ncol(ind_mat))[-1]) {
                ind1 = ind_mat[ir, ic-1] # previous column
                ind2 = ind_mat[ir, ic]   # current column
                v1 = small_mat[i[ind1], j[ind1]]
                v2 = small_mat[i[ind2], j[ind2]]
                if(v1 * v2 > 0) { # if they have the same sign
                    col = ifelse(v1 > 0, "violet", "blue")
                    grid.segments(x[ind1], y[ind1], x[ind2], y[ind2],
                        gp = gpar(col = col, lwd = 2))
                    grid.points(x[c(ind1, ind2)], y[c(ind1, ind2)], 
                        pch = 16, gp = gpar(col = col), size = unit(4, "mm"))
                }
            }
        }
        if(slice_r != slice_c) {
            grid.rect(gp = gpar(lwd = 2, fill = "transparent"))
        }
    }
)

p2 = grid.grabExpr(draw(ht, padding = unit(c(5, 5, 5, 5), "mm")))



library(cola)
rl = readRDS(url("https://jokergoo.github.io/cola_examples/TCGA_GBM/TCGA_GBM_subgroup.rds"))
res = rl["ATC:skmeans"]


m = get_matrix(rl)
sig = get_signatures(res, k = 4, plot = FALSE)
cl = get_classes(res, k = 4)[, 1]

m2 = m[sig$which_row, ]
m2 = t(scale(t(m2)))

col_fun = colorRamp2(c(-2, 0, 2), c("green", "white", "red"))

anno = get_anno(rl)
anno_col = get_anno_col(rl)

cl_col = cola:::brewer_pal_set2_col[1:4]
names(cl_col) = 1:4

ht = Heatmap(m2, name = "mat", show_row_names = FALSE, show_row_dend = FALSE, show_column_names = FALSE, 
    column_title = "C) Cluster all columns together",
    bottom_annotation = HeatmapAnnotation(df = anno, col = anno_col),
    top_annotation = HeatmapAnnotation(class = as.character(cl), col = list(class = cl_col)))
p3 = grid.grabExpr(draw(ht, merge_legend = TRUE))


ht = Heatmap(m2, name = "mat", row_split = sig$km, column_split = cl,
    column_title = "D) Columns are pre-split by 'class'",
    show_row_names = FALSE, show_row_dend = FALSE, show_column_names = FALSE, 
    bottom_annotation = HeatmapAnnotation(df = anno, col = anno_col),
    top_annotation = HeatmapAnnotation(class = as.character(cl), col = list(class = cl_col)))

p4 = grid.grabExpr(draw(ht, merge_legend = TRUE))

p_empty = rectGrob(gp = gpar(col = NA, fill = NA))

pdf("figure2.pdf", width = 10, height = 10)
library(cowplot)
p = plot_grid(plot_grid(p1, p2, nrow = 1),
              plot_grid(p_empty, p3, p4, nrow = 1, rel_widths = c(0.1, 1, 1)), 
              rel_heights = c(1.1, 1),
              nrow = 2)
print(p)
dev.off()



pdf("supplementary_figure1.pdf")
dimension_reduction(m2, col = cl_col[cl], method = "t-SNE")
dev.off()


