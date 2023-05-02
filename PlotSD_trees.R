source("plot_expr.r")

plot.expr.1 = function() {
#  png("image/expr/pha-4.png", width=1500, height=600)
  pdf("image/pha-4.pdf", width=11, height=6)

  col = rgb(scale.to.unit(expr.cell[246,]), 0, 0)
#    scale.to.unit(0 * expr.cell[246,]),
#    scale.to.unit(0))
  names(col) = names(expr.cell[246,])
  plot.segments.per.cell(col, cell.time.on.off)
  title("pha-4", cex.main=2)
  dev.off()
}
