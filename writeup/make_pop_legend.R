# Make an inset legend for the "corner pca" plots with space.
library(colorspace)

grid.width <- 4
grid.mat <- matrix(1:(grid.width^2), nrow=grid.width)
pop.cols <- rainbow_hcl(grid.width^2)
pop.pch <- seq_len(grid.width^2)

pdf(file="inset_pop_key.pdf", width=1, height=1, pointsize=10)
par(mar=c(0,0,0,0))
    plot(row(grid.mat), col(grid.mat),
         cex=3, pch=pop.pch, col=pop.cols, lwd=2,
         xlab='', ylab='', xaxt='n', yaxt='n',
         xlim=c(0.5, grid.width+0.5),
         ylim=c(0.5, grid.width+0.5) )
dev.off()
