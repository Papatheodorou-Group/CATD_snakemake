# Cellmix Installation script
#
# @zgr2788

suppressMessages(library(devtools))
suppressMessages(library(BiocManager))
devtools::install_github('shenorrlab/csSAM', auth_token = "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")

if (!require("GSEABase", quietly = TRUE)){
    BiocManager::install("GSEABase")
}

if (!require("Biobase", quietly = TRUE)){
    BiocManager::install("Biobase")
}

if (!require("genefilter", quietly = TRUE)){
    BiocManager::install("genefilter")
}

if (!require("preprocessCore", quietly = TRUE)){
    BiocManager::install("preprocessCore")
}


if (!require("NMF", quietly = TRUE)){
    install.packages("NMF", repos="http://cran.us.r-project.org")
}


if (!require("beeswarm", quietly = TRUE)){
    install.packages("beeswarm", repos="http://cran.us.r-project.org")
}

if (!require("quadprog", quietly = TRUE)){
    install.packages("quadprog", repos="http://cran.us.r-project.org")
}


if (!require("corpcor", quietly = TRUE)){
    install.packages("corpcor", repos="http://cran.us.r-project.org")
}

if (!require("gtools", quietly = TRUE)){
    install.packages("gtools", repos="http://cran.us.r-project.org")
}

if (!require("limSolve", quietly = TRUE)){
    install.packages("limSolve", repos="http://cran.us.r-project.org")
}

if (!require("matrixStats", quietly = TRUE)){
    install.packages("matrixStats", repos="http://cran.us.r-project.org")
}

if (!require("bibtex", quietly = TRUE)){
    install.packages("bibtex", repos="http://cran.us.r-project.org")
}

if (!require("CellMix", quietly = TRUE)){
    install.packages("Modules/CellMix/BiocInstaller_1.8.3.tar.gz", repos=NULL, type = "source")
    install.packages("Modules/CellMix/CellMix_1.6.2_R_x86_64-conda-linux-gnu.tar.gz", repos=NULL, type = "source")
}

suppressMessages(library(CellMix))
print("Installed CellMix successfully!")
