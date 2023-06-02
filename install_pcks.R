if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if(!require("cmapR", quietly = T)) {BiocManager::install("cmapR")}


if(!require("shinydashboard", quietly = T)) {install.packages("shinydashboard")}
if(!require("tidyverse", quietly = T)) {install.packages("tidyverse")}
if(!require("cowplot", quietly = T)) {install.packages("cowplot")}
if(!require("UpSetR", quietly = T)) {install.packages("UpSetR")}
if(!require("plotly", quietly = T)) {install.packages("plotly")}
if(!require("shinyWidgets", quietly = T)) {install.packages("shinyWidgets")}
if(!require("ggrepel", quietly = T)) {install.packages("ggrepel")}
if(!require("heatmaply", quietly = T)) {install.packages("heatmaply")}
if(!require("shinyShortcut", quietly = T)) {install.packages("shinyShortcut")}
if(!require("testthat", quietly = T)) {install.packages("testthat")}
if(!require("colourvalues", quietly = T)) {install.packages("colourvalues")}
if(!require("shinycssloaders", quietly = T)) {install.packages("shinycssloaders")}
if(!require("GGally", quietly = T)) {install.packages("GGally")}
if(!require("shinyjs", quietly = T)) {install.packages("shinyjs")}
if(!require("FactoMineR", quietly = T)) {install.packages("FactoMineR")}




# library(shinyShortcut)
# shinyShortcut(shinyDirectory = getwd(), OS = .Platform$OS.type, gitIgnore = FALSE)
