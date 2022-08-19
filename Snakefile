## Driver script for the deconvolution assessment pipeline
##
## @zgr2788
##
##
## Description:
## This file is essentially a driver for the different modules that are
## a part of the deconvolution assessment network.
##
##
##
##
##
##
##
##
##
##
##
##


#Helper functions to fetch inputs from other modules
def getVioPlots():
        metricsList = config['resMetrics']
        plotsList = [str("Plots/Vioplot_" + metric + ".png") for metric in metricsList if (metric == "rmse" or metric == "avgcos")]

        return(plotsList)

def getPlots():
        metricsList = config['resMetrics']
        plotsList = [str("Plots/Barplot_" + metric + ".png") for metric in metricsList]

        return(plotsList)



def getSelectedMetrics():
        metricsList = config['resMetrics']
        outputList = [str("Metrics/res_" + metric + ".rds") for metric in metricsList]

        return(outputList)



def getMethods():
        Methods = [str("Output/" + config['sampleName'] + "_res_" + methods + ".rds") for methods in config['deconMethods']]

        return(Methods)



def getBulks():
        inList = []

        if config['seededRun']:
            seedStatus = "_seeded"
        else:
            seedStatus = ""

        if config["stParam"]['scaleFirst']:
            st = "_scaled_transformed"
        else:
            st = "_transformed_scaled"

        filename_T = str("Input/Normalized_tables/" + config['sampleName'] + "_pbulks" + seedStatus + st + ".rds")
        filename_P = str("Input/Psuedobulks/" + config['sampleName'] + "_props" + seedStatus + ".rds")

        inList.append(filename_T)
        inList.append(filename_P)

        return inList



def getC2(inList):
        if config['seededRun']:
            seedStatus = "_seeded"
        else:
            seedStatus = ""

        filename_C2 = str("Input/References/" + config["sampleName"] + "_C2" + seedStatus + ".rds")

        inList.append(filename_C2)

        return inList



def getC1(inList):
        if config['seededRun']:
            seedStatus = "_seeded"
        else:
            seedStatus = ""

        filename_C1 = str("Input/References/" + config["sampleName"] + "_C1" + seedStatus + ".rds")

        inList.append(filename_C1)

        return inList



def getRefVar(inList):
        if config['seededRun']:
            seedStatus = "_seeded"
        else:
            seedStatus = ""

        filename_refVar = str("Input/References/" + config["sampleName"] + "_refVar" + seedStatus + ".rds")

        inList.append(filename_refVar)

        return inList



def getC0(inList):
        if config['seededRun']:
            seedStatus = "_seeded"
        else:
            seedStatus = ""

        if config["stParam"]['scaleFirst']:
            st = "_scaled_transformed"
        else:
            st = "_transformed_scaled"

        filename_C0 = str("Input/Normalized_tables/" + config["sampleName"] + "_C0" + seedStatus + st + ".rds")

        inList.append(filename_C0)

        return inList



def getPhenData(inList):
        if config['seededRun']:
            seedStatus = "_seeded"
        else:
            seedStatus = ""

        filename_phenData = str("Input/References/" + config["sampleName"] + "_phenData" + seedStatus + ".rds")

        inList.append(filename_phenData)

        return inList


configfile: 'config.yaml'
Methods = [str("Output/" + config['sampleName'] + "_res_" + methods + ".rds") for methods in config['deconMethods']]


#Prep metamodule
include: "Modules/Convert_split/Snakefile"

if not config['realBulk']:
    include: "Modules/Psuedobulk/Snakefile"
else:
    print("You chose to use real bulks.\nIf you see an error after this message, you either:\n-Forgot to place the pbulks and proportions\n-Placed them in the wrong place (should be under Input/Psuedobulks, typo intended)\n-The naming is wrong")

include: "Modules/C_generation/Snakefile"
include: "Modules/Scale_transform/Snakefile"

#Deconv metamodule
include: "Modules/debCAM/Snakefile"
include: "Modules/CDSeq/Snakefile"
include: "Modules/DeconRNASeq/Snakefile"
include: "Modules/OLS/Snakefile"
include: "Modules/NNLS/Snakefile"
include: "Modules/CIBERSORT/Snakefile"
include: "Modules/RLR/Snakefile"
include: "Modules/FARDEEP/Snakefile"
include: "Modules/DCQ/Snakefile"
include: "Modules/elasticNET/Snakefile"
include: "Modules/lasso/Snakefile"
include: "Modules/ridge/Snakefile"
include: "Modules/EPIC/Snakefile"
include: "Modules/dtangle/Snakefile"
include: "Modules/CellMix/Snakefile"
include: "Modules/ADAPTS/Snakefile"
include: "Modules/EpiDISH/Snakefile"
include: "Modules/MuSiC/Snakefile"
include: "Modules/BisqueRNA/Snakefile"
include: "Modules/SCDC/Snakefile"
include: "Modules/DWLS/Snakefile"
include: "Modules/bseqsc/Snakefile"
include: "Modules/CPM/Snakefile"
include: "Modules/TIMER/Snakefile"

#Results metamodule
include: "Modules/Res_explore/Snakefile"

#Grab all methods from config

#Target methods
rule all:
    input:
        getPlots(),
        getVioPlots(),
        "Benchmarks_summarized.png"
