#convertcolor.R

#takes a vector of a continuous variable and splits it into ncolor colors in a gradient from start to end

convertcolor <- function(v, ncolors = 10, percentile=FALSE, start="white",end="blue") {
  if (!require("RColorBrewer")) {
    install.packages("RColorBrewer")
    library(RColorBrewer)
  }
  scaled <- ifelse( rep(percentile,length(v)),
                    rank(v)/length(v),
                    (v - min(v)) / (max(v)-min(v))
  )
  colfunc <- colorRampPalette(c(start, end))
  binned <- cut(scaled, seq(0,1, by = 1/(ncolors)), include.lowest=TRUE, labels = colfunc(ncolors) )
  return (as.character(binned))
}
