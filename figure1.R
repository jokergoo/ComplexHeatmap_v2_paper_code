

library(grid)

single_heatmap_layout = function() {
th = convertHeight(grobHeight(textGrob("a")) * 3, "mm")
pushViewport(viewport(layout = grid.layout(nr = 3, nc = 3,
	width = unit.c(th*4, unit(1, "null"), th*4),
	height = unit.c(th*4, unit(1, "null"), th*4)),
	width = unit(1, "npc") - unit(4, "mm"),
	height = unit(1, "npc") - unit(4, "mm")))
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
name = c("column annotations", "column labels", "dendrogram", "title")
for(i in 1:4) {
	pushViewport(viewport(y = i/4, height = 1/4, just = "top"))
	grid.rect()
	grid.text(name[i])
	popViewport()
}
popViewport()

pushViewport(viewport(layout.pos.row = 3, layout.pos.col = 2))
name = rev(c("column annotations", "column labels", "dendrogram", "title"))
for(i in 1:4) {
	pushViewport(viewport(y = i/4, height = 1/4, just = "top"))
	grid.rect()
	grid.text(name[i])
	popViewport()
}
popViewport()

pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
name = rev(c("row annotations", "row labels", "dendrogram", "title"))
for(i in 1:4) {
	pushViewport(viewport(x = i/4, width = 1/4, just = "right"))
	grid.rect()
	grid.text(name[i], rot = 90)
	popViewport()
}
popViewport()

pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 3))
name = c("row annotations", "row labels", "dendrogram", "title")
for(i in 1:4) {
	pushViewport(viewport(x = i/4, width = 1/4, just = "right"))
	grid.rect()
	grid.text(name[i], rot = 90)
	popViewport()
}
popViewport()

pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 2))
name = matrix(c("matrix, row slice 2\ncolumn slice1", 
	         	"matrix, row slice 2\ncolumn slice2",
	         	"matrix, row slice 1\ncolumn slice1",
	         	"matrix, row slice 1\ncolumn slice2"), nrow = 2)
for(i in 1:2) {
	for(j in 1:2) {
		pushViewport(viewport(x = i/2, y = j/2, width = 1/2, height = 1/2, just = c("right", "top")))
		pushViewport(viewport(width = unit(1, "npc") - unit(4, "mm"),
			height = unit(1, "npc") - unit(4, "mm")))
		grid.rect(gp = gpar(col = "red", fill = "#FF000020"))
		grid.text(name[i, j])
		popViewport()
		popViewport()
	}
}
popViewport()

popViewport()

}


heatmap_annotation_layout = function() {

th = convertHeight(grobHeight(textGrob("a")) * 3, "mm")
pushViewport(viewport(
	width = unit(1, "npc") - unit(20, "mm"),
	height = unit(1, "npc") - unit(4, "mm")))
grid.rect()
for(i in 1:3) {
	name = c("sinlge annotation 1", "single annotation 2", "single annotation 3")

	pushViewport(viewport(x = i/3, width = 1/3, just = "right"))
	pushViewport(viewport(width = unit(1, "npc") - unit(4, "mm"),
		height = unit(1, "npc") - unit(4, "mm")))
	grid.rect(gp = gpar(col = "blue", fill = "#0000FF20"))
	grid.text(name[i], rot = 90)
	popViewport()
	popViewport()
}
popViewport()

}



heatmap_list_layout = function(direction = "horizontal") {

th = convertHeight(grobHeight(textGrob("a")) * 3, "mm")
pushViewport(viewport(layout = grid.layout(nr = 3, nc = 3,
	width = unit.c(th*2, unit(1, "null"), th*2),
	height = unit.c(th*2, unit(1, "null"), th*2)),
	width = unit(1, "npc") - unit(4, "mm"),
	height = unit(1, "npc") - unit(4, "mm")))
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
name = c("title", "legends")
for(i in 1:2) {
	pushViewport(viewport(y = i/2, height = 1/2, just = "top"))
	grid.rect()
	grid.text(name[i])
	popViewport()
}
popViewport()

pushViewport(viewport(layout.pos.row = 3, layout.pos.col = 2))
name = rev(c("title", "legends"))
for(i in 1:2) {
	pushViewport(viewport(y = i/2, height = 1/2, just = "top"))
	grid.rect()
	grid.text(name[i])
	popViewport()
}
popViewport()

pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
name = rev(c("title", "legends"))
for(i in 1:2) {
	pushViewport(viewport(x = i/2, width = 1/2, just = "right"))
	grid.rect()
	grid.text(name[i], rot = 90)
	popViewport()
}
popViewport()

pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 3))
name = c("title", "legends")
for(i in 1:2) {
	pushViewport(viewport(x = i/2, width = 1/2, just = "right"))
	grid.rect()
	grid.text(name[i], rot = 90)
	popViewport()
}
popViewport()

pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 2))
name = c("heatmap1", "heatmap2", "row annotation")
for(i in 1:3) {
	if(direction == "horizontal") {
		name = c("heatmap 1", "heatmap 2", "row annotation")

		pushViewport(viewport(x = i/3, width = 1/3, just = "right"))
		pushViewport(viewport(width = unit(1, "npc") - unit(4, "mm"),
			height = unit(1, "npc") - unit(4, "mm")))
		grid.rect(gp = gpar(col = "green", fill = "#00FF0020"))
		grid.text(name[i], rot = 90)
		popViewport()
		popViewport()
	} else {
		name = rev(c("heatmap 1", "heatmap 2", "column annotation"))

		pushViewport(viewport(y = i/3, height = 1/3, just = "top"))
		pushViewport(viewport(width = unit(1, "npc") - unit(4, "mm"),
			height = unit(1, "npc") - unit(4, "mm")))
		grid.rect(gp = gpar(col = "red", fill = "#FF000020"))
		grid.text(name[i])
		popViewport()
		popViewport()
	}
}
popViewport()

popViewport()

}


pdf("figure1.pdf", width = 11.3, height = 4.4)

grid.newpage()

th = convertHeight(grobHeight(textGrob("a")) * 3, "mm")

pushViewport(viewport(layout = grid.layout(nr = 2, nc = 3, 
	height = unit.c(th, unit(1, "null")), 
	widths = unit(c(2.5, 1, 1.2), "null")),
	gp = gpar(fontsize = 10)))
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
grid.text("A) A single heatmap", gp = gpar(fontsize = 12))
popViewport()

pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
single_heatmap_layout()
popViewport()

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
grid.text("B) Heatmap annotation", gp = gpar(fontsize = 12))
popViewport()

pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 2))
heatmap_annotation_layout()
popViewport()

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 3))
grid.text("C) Heatmap list", gp = gpar(fontsize = 12))
popViewport()

pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 3))
heatmap_list_layout()
popViewport()

popViewport()

dev.off()


