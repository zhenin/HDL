load("Figure3.rda")

### plot category for all proteins, panel a
a <- nrow(new) ## NEW, num of protein-trait pairs with proteins as new target
b <- sum(rr$Category == 0) ## Druggable,
c <- sum(rr$Category == 2) ## re-purposing
e <- sum(rr$Category == 1) ## matched+
f <- sum(rr$Category == -1) ## matched-
pa <- c(a, b ,c, e, f)



nejm.colors <- function (n = 8, alpha = 1){
  pal <- c("#BC3C29", "#0072B5", "#E18727", "#20854E", "#7876B1",
           "#6F99AD", "#FFDC91", "#EE4C97")
  acode <- substr(rgb(0, 0, 0, alpha = alpha), 8, 9)
  return(paste0(pal, acode)[1:n])
}

colora <- nejm.colors(8, .9)[1:4]


### category by protein
colorb <- colora[-1]
bb <- t(table(rr$protein, rr$Category))

bb <- bb[c('0', '2', '1'),]

pdf('Figure3.pdf', width = 16*.9, height = 9*0.9)
par(mar = c(7,4,2,0), xpd = TRUE, mgp = c(2.2,.8,0))
layout(matrix(c(1,2,3,3), ncol = 2, byrow = TRUE), widths = c(1,5), heights = c(4,5))

#a
bar1 <- barplot(pa, col = colora, ylim = c(0, max(pa)+35*0.5), xaxt = 'n', ylab = '#Pairs', cex.axis = 1.2, cex.lab = 1.2)
text(bar1, pa+max(pa)*0.08, paste(pa, sep=""), cex = 1)
text(bar1, -max(pa)*0.08, c('New', 'Druggable', 'Re-purposing', expression('Matched'^'+'), expression('Matched'^'-')), srt = 60, adj = 1, cex = 1.2)
#b
par(mar = c(7, 7, 2, 0))
bar2 <- barplot((bb), beside = TRUE, width = 10, cex.lab = 1.2,
                col = colorb, ylim = c(0, max(bb)*1.2), cex.axis = 1.2, #plot = FALSE,
                xaxt = 'n', ylab = '#Pairs')

text(bar2, (bb)+max(bb)*.05, paste((bb), sep=""), cex = 1)
text(colMeans(bar2), -max(bb)*0.08, colnames(bb), srt = 60, adj = 1, cex = 1.2)

dev.off()