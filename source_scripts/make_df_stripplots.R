# function that takes as input the Affymetrix probeset ID(s) and gives as output a data.frame
prepareData <- function(x){
  selected_genes <- t(prDat[x, ])
  pDat <- cbind(selected_genes, prDes)
  pDat <- with(pDat, data.frame(sidChar, sidNum, devStage, gType, 
          gene = factor(rep(c(colnames(pDat[x])), each = nrow(pDat))), 
          gExp = c(selected_genes)))
}

# function to plot stripplot using lattice
makeStripplot <- function (x, ...) {
  stripplot(gExp ~ devStage | gene, x,
            group = gType, jitter.data = TRUE,
            auto.key = TRUE, type = c('p', 'a'), grid = TRUE, ...)
}

# function to plot stripplot using ggplot2
makeStripplotGg <- function (x) {
  ggplot(x, aes(x = devStage, y = gExp, color = gType, group = gType)) +
    geom_point() + stat_smooth(se = F) + facet_wrap(~ gene)
}