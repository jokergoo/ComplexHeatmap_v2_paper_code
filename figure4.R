library(ComplexHeatmap)
library(circlize)

res_list = readRDS("~/project/development/ComplexHeatmap-reference/data/meth.rds")
type = res_list$type
mat_meth = res_list$mat_meth
mat_expr = res_list$mat_expr
direction = res_list$direction
cor_pvalue = res_list$cor_pvalue
gene_type = res_list$gene_type
anno_gene = res_list$anno_gene
anno_enhancer = res_list$anno_enhancer

dist = sapply(1:nrow(mat_meth), function(i) {
	if(mean(mat_meth[i, ]) < 0.2) {
		runif(1, 0, 1000)
	} else if(mean(mat_meth[i, ]) < 0.4) {
		runif(1, 1000, 5000)
	} else {
		runif(1, 1000, 50000)
	}
})

column_tree = hclust(dist(t(mat_meth)))
column_order = column_tree$order

library(RColorBrewer)
meth_col_fun = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))
direction_col = c("hyper" = "red", "hypo" = "blue")
expr_col_fun = colorRamp2(c(-2, 0, 2), c("green", "white", "red"))
pvalue_col_fun = colorRamp2(c(0, 2, 4), c("white", "white", "red"))
gene_type_col = structure(brewer.pal(length(unique(gene_type)), "Set3"), 
    names = unique(gene_type))
anno_gene_col = structure(brewer.pal(length(unique(anno_gene)), "Set1"), 
    names = unique(anno_gene))
dist_col_fun = colorRamp2(c(0, 10000), c("black", "white"))
enhancer_col_fun = colorRamp2(c(0, 1), c("white", "orange"))

ha = HeatmapAnnotation(type = type, 
    col = list(type = c("Tumor" = "pink", "Control" = "royalblue")),
    annotation_name_side = "left", 
    annotation_legend_param = list(type = list(direction = "horizontal")))
ha2 = HeatmapAnnotation(type = type, 
    col = list(type = c("Tumor" = "pink", "Control" = "royalblue")), 
    show_legend = FALSE)

ht_list = Heatmap(mat_meth, name = "methylation", col = meth_col_fun,
    column_order= column_order,
    top_annotation = ha, column_title = "Methylation",
    heatmap_legend_param = list(direction = "horizontal", legend_width = unit(3, "cm"))) +
    Heatmap(direction, name = "direction", col = direction_col,
    	heatmap_legend_param = list(direction = "horizontal")) +
    Heatmap(mat_expr[, column_tree$order], name = "expression", 
        col = expr_col_fun, 
        column_order = column_order, 
        top_annotation = ha2, column_title = "Expression",
        heatmap_legend_param = list(direction = "horizontal", legend_width = unit(3, "cm"))) +
    Heatmap(cor_pvalue, name = "p-value", col = pvalue_col_fun,
    	heatmap_legend_param = list(direction = "horizontal", at = c(0, 2, 4), labels = c("1", "0.01", "<1e-4"))) +
    Heatmap(gene_type, name = "gene type", col = gene_type_col,
    	heatmap_legend_param = list(direction = "horizontal",  nrow = 2)) +
    Heatmap(anno_gene, name = "gene annotation", col = anno_gene_col,
    	heatmap_legend_param = list(direction = "horizontal", nrow = 2)) +
    rowAnnotation(dist_to_TSS = anno_points((dist+1), width = unit(4, "cm"))) +
    Heatmap(anno_enhancer, name = "Enhancer", col = enhancer_col_fun, 
        cluster_columns = FALSE, column_title = "Enhancer",
        heatmap_legend_param = list(direction = "horizontal", legend_width = unit(3, "cm")))

p1 = grid.grabExpr(draw(ht_list, row_km = 2, row_split = direction,
    column_title = "A) Comprehensive correspondence between methylation, expression and other genomic features", 
    merge_legends = TRUE, heatmap_legend_side = "bottom"), width = 12)



lt = readRDS(system.file("extdata", package = "ComplexHeatmap", "dmr_summary.rds"))
attach(lt)

library(circlize)
meth_col_fun = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))
corr_col = c("green", "red")
dist_tss_col = c("#FF0000", "#FF7352", "#FFB299", "#FFD9CB")
gene_anno_col = c("green", "blue")
cgi_anno_col = c("#FFA500", "#FFD191")
z_score_col_fun = colorRamp2(c(-200, 0, 200), c("green", "white", "red"))
state_col = c("#FF0000", "#008000", "#C2E105", "#8A91D0", "#CD5C5C", "#808080", "#000000")

anno_width = unit(3.4, "cm")
ht_list = rowAnnotation(text = anno_text(label, location = unit(1, "npc"), just = "right", 
    gp = gpar(fontsize = 12)))

ht_list = ht_list + Heatmap(mean_meth, name = "mean_meth", col = meth_col_fun, 
    cluster_rows = FALSE, row_title = NULL, cluster_columns = FALSE, show_row_names = FALSE, column_names_rot = 45,
    column_names_gp = gpar(fontsize = 9),
    heatmap_legend_param = list(title = "Methylation", direction = "horizontal", legend_width = unit(3, "cm")), 
    width = ncol(mean_meth)*unit(4, "mm")) +
rowAnnotation("n_gr" = anno_barplot(n_gr, bar_width = 1, width = anno_width), 
    show_annotation_name = FALSE) +
rowAnnotation("n_corr" = anno_barplot(n_corr, bar_width = 1, gp = gpar(fill = corr_col), 
    width = anno_width), show_annotation_name = FALSE) +
rowAnnotation("dist_tss" = anno_barplot(dist_tss, bar_width = 1, gp = gpar(fill = dist_tss_col), 
    width = anno_width), show_annotation_name = FALSE) +
rowAnnotation("gene_anno" = anno_barplot(gene_anno, bar_width = 1, gp = gpar(fill = gene_anno_col), 
    width = anno_width), show_annotation_name = FALSE) +
rowAnnotation("cgi_anno" = anno_barplot(cgi_anno, bar_width = 1, gp = gpar(fill = cgi_anno_col), 
    width = anno_width), show_annotation_name = FALSE) +
Heatmap(mat_enrich_gf, name = "enrich_gf", col = z_score_col_fun, cluster_columns = FALSE,
    width = unit(ncol(mat_enrich_gf)*4, "mm"), column_title = "", column_names_rot = 45,
    column_names_gp = gpar(fontsize = 9),
    heatmap_legend_param = list(title = "Z-score",direction = "horizontal", legend_width = unit(4, "cm"))) +
rowAnnotation("pct_st" = anno_barplot(mat_pct_st, bar_width = 1, gp = gpar(fill = state_col), 
    width = anno_width), show_annotation_name = FALSE) +
Heatmap(mat_enrich_st, name = "enrich_st", col = z_score_col_fun, cluster_columns = FALSE, 
    width = unit(ncol(mat_enrich_st)*4, "mm"), column_title = "", show_heatmap_legend = FALSE,
    column_names_gp = gpar(col = state_col, fontsize = 9), show_row_names = FALSE, column_names_rot = 45)

lgd_list = list(
    Legend(labels = c("gene", "intergenic"), title = "Gene annotation", 
        legend_gp = gpar(fill = gene_anno_col)),
    Legend(labels = c("<1kb", "1kb~5kb", "5kb~10kb", ">10kb"), title = "Distance to TSS", 
        legend_gp = gpar(fill = dist_tss_col), nrow = 2),
    Legend(labels = c("CGI", "CGI shore"), title = "CGI annotation", 
        legend_gp = gpar(fill = cgi_anno_col)),
    Legend(labels = colnames(mat_enrich_st), title = "Chromatin states", 
        legend_gp = gpar(fill = state_col), nrow = 2)
)


p2 = grid.grabExpr({
	draw(ht_list, padding = unit(c(2, 2, 20, 2), "mm"), row_split = gsub("\\d+$", "", label), 
	    heatmap_legend_list = lgd_list, heatmap_legend_side = "bottom")
	anno_title = c("n_gr" = "Number of\nDMRs", "n_corr" = "Significantly\ncorrelated genes",
	    "gene_anno" = "Gene annotation", "dist_tss" = "Distance to TSS",
	    "cgi_anno" = "CGI annotation", "pct_st" = "Overlap to\nChromatin states")
	for(an in names(anno_title)) {
	    decorate_annotation(an, {
	        grid.text(anno_title[an], y = unit(1, "npc") + unit(3, "mm"), just = "bottom")
	    })
	}
	ht_title = c("mean_meth" = "Mean\nmethylation", "enrich_gf" = "Enrichment to\ngenomic features",
	    "enrich_st" = "Enrichment to\nchromatin states")
	for(an in names(ht_title)) {
	    decorate_heatmap_body(an, {
	        grid.text(ht_title[an], y = unit(1, "npc") + unit(3, "mm"), just = "bottom")
	    })
	}
}, width = 12)

library(cowplot)

p_empty = rectGrob(gp = gpar(col = NA, fill = NA))
pdf("figure4.pdf", width = 13, height = 10)
p = plot_grid(
    plot_grid(p_empty, p1, p_empty, nrow = 1, rel_widths = c(0.15, 1, 0.15)), 
    p2, nrow = 2, rel_heights = c(1.5, 1))
print(p)
dev.off()

