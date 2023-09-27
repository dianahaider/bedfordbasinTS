install.packages("UpSetR")
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("ComplexHeatmap")

library(UpSetR)
library(ComplexHeatmap)


#upload data
upsetplots <- read_csv("data/upsetplots.csv", 
                       +     col_types = cols(Phylum = col_skip()))

set.seed(123)
m1 = make_comb_mat(upsetplots)
m1


upset(movies,attribute.plots=list(gridrows=60,plots=list(list(plot=scatter_plot, x="ReleaseDate", y="AvgRating"),
                                                         list(plot=scatter_plot, x="ReleaseDate", y="Watches"),list(plot=scatter_plot, x="Watches", y="AvgRating"),
                                                         list(plot=histogram, x="ReleaseDate")), ncols = 2))