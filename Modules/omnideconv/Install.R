# Omnideconv installation script
#
# @zgr2788

#if (!require("pak", quietly = TRUE)){
#  install.packages("pak", repos="http://cran.us.r-project.org")
#}

# install the `pak` package manager
install.packages("pak",repos = "http://cran.us.r-project.org")

# minimal installation
#suppressMessages(library(pak))
#pak::pkg_install("omnideconv/omnideconv")

# complete installation, including Python dependencies
suppressMessages(library(pak))
pak::pkg_install("omnideconv/omnideconv")

suppressMessages(library(omnideconv))
install.packages("e1071", repos="http://R-Forge.R-project.org")
suppressMessages(library(e1071))

print("Installed omnideconv successfully!")

#Manual
#mamba install r-base geos r-gert r-vctrs "gdal=3.3.3.""
