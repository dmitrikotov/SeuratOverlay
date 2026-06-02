# SeuratOverlay
R package to overlay cluster metadata on top of the Seurat FeaturePlot function

# Installing SeuratOverlay
To install SeuratOverlay and the required depndencies from github run the following code:

if(!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools") 
}

devtools::install_github("dmitrikotov/SeuratOverlay")
