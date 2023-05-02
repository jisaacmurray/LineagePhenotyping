# Various figures.

library(grid)

od = setwd("/murrlab/jburdick/gcb/work/svnroot/trunk/deconvolve/R/")
# source("import.r")
# source("tree/plot_expr.r")
source("deconvolve.r")
setwd(od)

source("R/plot/plot_expr.r")

# Some expression graphs.
plot.expr.1 = function() {
#  png("image/expr/pha-4.png", width=1500, height=600)
  pdf("image/expr/pha-4.pdf", width=11, height=6)

  col = rgb(scale.to.unit(expr.cell[246,]), 0, 0)
#    scale.to.unit(0 * expr.cell[246,]),
#    scale.to.unit(0))
  names(col) = names(expr.cell[246,])
  plot.segments.per.cell(col, cell.time.on.off)
  title("pha-4", cex.main=2)
  dev.off()
}


# Plots of distinguishable regions.
show.regions = function() {




}

# Plots the cell-cell correlation matrix, including the lineage.
# FIXME not yet working.
plot.with.lineage = function(x) {
   grey.colors = hsv(0, 0, c(0:255)/255)

  pdf("R/plot/cellCellCorrelation.pdf", width=5, height=7.5)
  # ??? should I be using grid graphics?
  def.par = par(no.readonly=TRUE)
  layout(matrix(c(1,2), nrow=2), widths=c(1), heights=c(0.2, 0.8))

#  par(mar=c(0,4,1,0)
  col = rgb(scale.to.unit(0 * expr.cell[1,]), 0, 0)
#    scale.to.unit(0 * expr.cell[246,]),
#    scale.to.unit(0))
  names(col) = names(expr.cell[1,])
  plot.segments.per.cell(col)

  image(x, col=grey.colors, xaxt="n", yaxt="n", xlab="", ylab="", useRaster=TRUE)
# ??? use grid.raster(as.matrix(x), col=grey.colors)

if (FALSE) {

  loc.12.cell = sapply(lin.12.cell, function(x) cell.to.column[[x]])
  axis(1, at=loc.12.cell / 1341, label=lin.12.cell, cex.axis=2)
  axis(2, at=loc.12.cell / 1341, label=lin.12.cell, cex.axis=2)

  # add in rectangle around some lineages
  c.ABprp = cell.to.column[c(lin.list[["ABprp"]], recursive=TRUE)]
  c.ABplp = cell.to.column[c(lin.list[["ABplp"]], recursive=TRUE)]
  rect(min(c.ABprp)/1341, min(c.ABplp)/1341, max(c.ABprp)/1341, max(c.ABplp)/1341,
    border="blue", col=NULL, lwd=5)
  rect(min(c.ABplp)/1341, min(c.ABprp)/1341, max(c.ABplp)/1341, max(c.ABprp)/1341,
    border="blue", col=NULL, lwd=5)
}

  par(def.par)
  dev.off()
}

# XXX hacky conversion
# cor.cell.shrunk.1 = matrix(as.vector(cor.cell.shrunk), ncol=1341)
# system.time(plot.with.lineage(cor.cell.shrunk.1))



