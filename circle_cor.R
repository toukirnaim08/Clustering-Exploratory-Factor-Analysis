circle_cor = function(cor, axes = FALSE, xlab = "", title = "Correlation Plot",
     ylab = "", asp = 1, cex.lab = par("cex.lab"), cex = 0.75 * par("cex"),
     ...) {

# Originally written by Taiyun Wei
# Modified by Belinda Chiera (2012)
#     	     
# This function plots correlations as circles.  The larger the cirle, the higher the correlation.
# Circles that are black correspond to correlations > 0.
# Circles that are white correspond to correlations <= 0.	
# An example of how to run this function is as follows:
# In R, at the command prompt, type:
#
# source("circle_cor.R")
#
# and press enter.  If no error message appears, you can run the function.  Assuming you have the
# correlations stored in a matrix called corr.ytube you would call circle_cor as follows:
#
# > circle_cor(corr.ytube)
#
# and press enter and a plot should appear!  
# If you would like to change the name of the title of this plot, at the very top of this file there is the argument
# title = "Corelation Plot"
#
# just change the text Correlation Plot to whatever you would like to call the plot, ensuring the text is inside " " marks.
#
# If you would like to change the colours of the circles, scroll down towards the bottom of this file.  Look for these lines of code:
#
#bg[cor > 0] = "black"
#bg[cor <= 0] = "white"  
#
# Replace the words "black" and "white" with whichever colours you like (ensuring you write the colour with " " marks).  Then at
# the R command prompt, run the source command again and then run circle_cor and that's it!

rowdim <- dim(cor)[1]
coldim <- dim(cor)[2]
rowlabs <- dimnames(cor)[[1]]
collabs <- dimnames(cor)[[2]]
cols <- 1:coldim
rows <- 1:rowdim

if (is.null(rowlabs))
    rowlabs <- 1:rowdim
if (is.null(collabs))
    collabs <- 1:coldim
rowlabs <- as.character(rowlabs)
collabs <- as.character(collabs)
maxdim <- max(length(rows), length(cols))
plt <- par("plt")
 
xlabwidth <- max(strwidth(rowlabs[rows], units = "figure", cex = cex.lab))/(plt[2] - plt[1])
xlabwidth <- xlabwidth * maxdim/(1 - xlabwidth)
ylabwidth <- max(strwidth(collabs[cols], units = "figure", cex = cex.lab))/(plt[4] - plt[3])
ylabwidth <- ylabwidth * maxdim/(1 - ylabwidth)
 	
 n = nrow(cor)
 par(mar = c(0, 0, 2, 0), bg = "white")
 plot(c(-xlabwidth - 0.5, maxdim + 0.5), c(0.5, maxdim + 1 + ylabwidth), axes = axes, xlab = "",
         ylab = "", asp = 1, type = "n")
 ##add grid
 segments(rep(0.5, n + 1), 0.5 + 0:n, rep(n + 0.5, n + 1),
         0.5 + 0:n, col = "gray")
 segments(0.5 + 0:n, rep(0.5, n + 1), 0.5 + 0:n, rep(n + 0.5,
                 n), col = "gray")
 ##define circles' background color.
 ##black for positive correlation coefficient and white for negative
 bg = cor
 bg[cor > 0] = "blue"
 bg[cor <= 0] = "green"   
 
  
# plot n*n circles using vector language, suggested by Yihui Xie  
symbols(rep(1:n, each = n), rep(n:1, n), add = TRUE, inches = F, circles = as.vector(sqrt(abs(cor))/2), bg = as.vector(bg))  
text(rep(0, length(rows)), length(rows):1, labels = rowlabs[rows], adj = 1, cex = cex.lab)
text(cols, rep(length(rows) + 1, length(cols)), labels = collabs[cols], srt = 90, adj = 0, cex = cex.lab)
mtext(xlab, 1, 0)
mtext(ylab, 2, 0)
title(title)
} 


