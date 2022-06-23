

library(GetoptLong)
library(ComplexHeatmap)

set.seed(123)
image_png = sample(dir("~/project/development/ComplexHeatmap-reference/IcoMoon-Free-master/PNG/64px", full.names = TRUE), 20)

plot1 = function() {
	n = 7
	m = matrix(rnorm(100*n), nc = 100)
	m2 = matrix(rnorm(100*n), nc = 100)

	anno_pct = function(x) {

	    max_x = max(x)
	    text = paste0(sprintf("%.2f", x*100), "%")
	    cell_fun_pct = function(i) {
	        pushViewport(viewport(xscale = c(0, max_x)))
	        grid.roundrect(x = unit(1, "npc"), width = unit(x[i], "native"), height = unit(1, "npc") - unit(4, "pt"), 
	            just = "right", gp = gpar(fill = "#0000FF80", col = NA))
	        grid.text(text[i], x = unit(1, "npc"), just = "right")
	        popViewport()
	    }

	    AnnotationFunction(
	        cell_fun = cell_fun_pct,
	        var_import = list(max_x, x, text), 
	        which = "row",
	        width = max_text_width(text)*1.5
	    )
	}


	ha1 = rowAnnotation(
		"anno_simple" = runif(n),
		"anno_image" = anno_image(image_png[1:n], gp = gpar(col = "black"), space= unit(2, "mm"), width = unit(1, "cm")),
		"anno_points" = anno_points(matrix(runif(2*n), nc = 2), pch = 16:17, gp = gpar(col = 2:3), width = unit(2, "cm")),
		"anno_lines" = anno_lines(cbind(runif(n), runif(n)), gp = gpar(col = 2:3, lty = 1:2), width = unit(2, "cm")),
		"anno_loess" = anno_lines(runif(n), smooth = TRUE, width = unit(2, "cm")),
		"anno_barplot" = anno_barplot(cbind(runif(n), runif(n)), gp = gpar(fill = 2:3, col = 2:3), width = unit(2, "cm")),
		"anno_pct" = anno_pct(runif(n)),
		"anno_boxplot" = anno_boxplot(m, gp = gpar(fill = 1:n), width = unit(3.5, "cm")),
		"anno_text" = anno_text(gt_render(qq("row<sub>@{1:n}</sub>", collapse = FALSE), align_widths = TRUE, halign = 0.5, 
	                        r = unit(4, "pt"),
	                        padding = unit(c(2, 8, 2, 8), "pt")), 
	                    gp = gpar(box_col = "blue", box_lwd = 2, col = 1:10, fontfamily = "Times"), 
	                    just = "right", 
	                    location = unit(1, "npc")
	    ),
	    "anno_histogram" = anno_histogram(m2, gp = gpar(fill = 1:n), border = TRUE, width = unit(3.5, "cm")),
		"anon_violin" = anno_density(m, type = "violin", gp = gpar(fill = 1:10, col = c("white", rep("black", 7), "white", "black")), width = unit(3.5, "cm")),
		show_legend = FALSE, simple_anno_size = unit(1, "cm"), gap = unit(2, "mm"), show_annotation_name = FALSE
	)


	m = matrix(rnorm(1500), nc = 15)
	lt1 = apply(m, 2, function(x) data.frame(density(x)[c("x", "y")]))

	lt2 = lapply(1:15, function(x) cumprod(1 + runif(1000, -x/100, x/100)) - 1)

	# 20 rows
	ha3 = rowAnnotation(
		"anno_joyplot" = anno_joyplot(lt1, width = unit(4, "cm"), gp = gpar(fill = 1:10), transparency = 0.75, scale = 2),
		"anno_horizon" = anno_horizon(lt2),
	    gap = unit(2, "mm"), show_annotation_name = FALSE
	)


	grid.newpage()
	grid.rect(gp = gpar(fill = "white", col = NA))
	pushViewport(viewport(width = ComplexHeatmap:::width(ha1) + ComplexHeatmap:::width(ha3)))

	pushViewport(viewport(x = 0, y = unit(1, "npc") - unit(5, "mm"), width = ComplexHeatmap:::width(ha1), height = unit(1, "npc") - unit(2, "cm"), just = c("left", "top")))
	draw(ha1)
	decorate_annotation("anno_image", {
		grid.text("Images", y = unit(1, "npc") + unit(6, "mm"), just = "bottom", gp = gpar(fontsize = 12))
	})	
	decorate_annotation("anno_points", {
		grid.text("Points", y = unit(1, "npc") + unit(6, "mm"), just = "bottom", gp = gpar(fontsize = 12))
	})	
	decorate_annotation("anno_lines", {
		grid.text("Lines", y = unit(1, "npc") + unit(6, "mm"), just = "bottom", gp = gpar(fontsize = 12))
	})	
	decorate_annotation("anno_loess", {
		grid.text("Smoothed\nlines", y = unit(1, "npc") + unit(3, "mm"), just = "bottom", gp = gpar(fontsize = 12))
	})	
	decorate_annotation("anno_barplot", {
		grid.text("Barplot", y = unit(1, "npc") + unit(6, "mm"), just = "bottom", gp = gpar(fontsize = 12))
	})	
	decorate_annotation("anno_pct", {
		grid.text("Percent\nwith bars", y = unit(1, "npc") + unit(3, "mm"), just = "bottom", gp = gpar(fontsize = 12))
	})	
	decorate_annotation("anno_boxplot", {
		grid.text("Boxplot", y = unit(1, "npc") + unit(6, "mm"), just = "bottom", gp = gpar(fontsize = 12))
	})	
	decorate_annotation("anno_text", {
		grid.text("Customized\ntext", y = unit(1, "npc") + unit(3, "mm"), just = "bottom", gp = gpar(fontsize = 12))
	})	
	decorate_annotation("anno_histogram", {
		grid.text("Histogram", y = unit(1, "npc") + unit(6, "mm"), just = "bottom", gp = gpar(fontsize = 12))
	})
	decorate_annotation("anon_violin", {
		grid.text("Violin plot", y = unit(1, "npc") + unit(6, "mm"), just = "bottom", gp = gpar(fontsize = 12))
	})
	popViewport()

	pushViewport(viewport(x = ComplexHeatmap:::width(ha1) + unit(2, "mm"), y = unit(1, "npc") - unit(5, "mm"), width = ComplexHeatmap:::width(ha3) - unit(2, "mm"), height = unit(1, "npc") - unit(2, "cm"), just = c("left", "top")))
	draw(ha3)
	decorate_annotation("anno_joyplot", {
		grid.text("Joy plot", y = unit(1, "npc") + unit(6, "mm"), just = "bottom", gp = gpar(fontsize = 12))
	})
	decorate_annotation("anno_horizon", {
		grid.text("Horizon chart", y = unit(1, "npc") + unit(6, "mm"), just = "bottom", gp = gpar(fontsize = 12))
	})
	popViewport()

	popViewport()
}


p1 = grid.grabExpr(plot1())

nr1 = 40; nr2 = 80; nr3 = 60; nr = nr1 + nr2 + nr3
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
mat = mat[sample(nr, nr), ] # random shuffle rows and columns

plot2 = function() {
	genes = sample(nr, 20)
	labels = paste0("gene", 1:20)

	ht = Heatmap(mat, name = "mat", right_annotation = rowAnnotation(mark = anno_mark(at = genes, labels = labels)))
	draw(ht, column_title = "B) Mark annotation")
}

p2 = grid.grabExpr(plot2(), height = 10*1.2/2.3, width = 15*2/8)

library(ggplot2)

plot3 = function() {
	km = kmeans(mat, centers = 3)$cluster
	ind_list = split(1:nrow(mat), km)
	groups = c(rep("group1", 6), rep("group2", 8), rep("group3", 10))
	ht = Heatmap(mat, name = "mat", row_split = km,
		column_split = groups, cluster_column_slices = FALSE,
		right_annotation = rowAnnotation(link = anno_link(align_to = ind_list,
			panel_fun = function(index) {
				mm = mat[index, ]
				x1 = rowMeans(mm[, 1:6])
				x2 = rowMeans(mm[, 6+1:8])
				x3 = rowMeans(mm[, 14+1:10])
				df = data.frame(value = c(x1, x2, x3), group = paste0("group", rep(1:3, each = nrow(mm))))
				p = ggplot(df, aes(group, value)) + geom_boxplot() + theme(axis.title.x=element_blank())
				p = grid.grabExpr(print(p))
				grid.draw(p)
				grid.segments(c(0, 1, 1), c(0, 0, 1), c(1, 1, 0), c(0, 1, 1))
			}, size = unit(3, "cm"), width = unit(7, "cm"), internal_line = FALSE, gap = unit(4, "mm")))
	)
	draw(ht, column_title = "C) Link annotation")
}

p3 = grid.grabExpr(plot3(), height = 10*1.2/2.3, width = 15*3/8)

plot4 = function() {
	library(circlize)
	mat = matrix(rnorm(100*10), nrow = 100)

	split = sample(letters[1:8], 100, replace = TRUE)
	text = lapply(unique(split), function(x) {
	    data.frame(paste0("functioin", 1:8), col = rand_color(8, friendly = TRUE), fontsize = runif(8, 6, 14))
	})
	names(text) = unique(split)

	ht = Heatmap(mat, name = "mat", cluster_rows = FALSE, row_split = split,
	    right_annotation = rowAnnotation(wc = anno_textbox(split, text, max_width = unit(8, "cm")))
	)
	draw(ht, column_title = "D) Textbox annotation")
}

p4 = grid.grabExpr(plot4(), height = 10*1.2/2.3, width = 15*3/8)

library(cowplot)

p_empty = textGrob("A) Different annotation graphics supported in ComplexHeatmap", gp = gpar(fontsize = 12))

pdf("figure3.pdf", width = 15, height = 10)
print(plot_grid(p_empty, p1, 
	plot_grid(p2, p3, p4, nrow = 1, rel_widths = c(1.5, 3, 3)), 
	nrow = 3, rel_heights = c(0.2, 0.9, 1.2)))
dev.off()

