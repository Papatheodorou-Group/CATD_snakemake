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

configfile: 'config.yaml'
include: "Modules/Convert_split/Snakefile"
include: "Modules/Psuedobulk/Snakefile"
include: "Modules/C_generation/Snakefile"


outFile = list()
inFile = config['Input']

if config['seededRun'] == False:
    pbulks = inFile + '_pbulks.rds'
    props = inFile + '_props.rds'

else:
    pbulks = inFile + '_pbulks_seeded.rds'
    props = inFile + '_props_seeded.rds'

outFile.append(pbulks)
outFile.append(props)

rule normalizeMatrices:
    input:
        outFile,
        'Input/References/Hrvatin_afteint_C1.rds',
        'Input/References/Hrvatin_afteint_C2.rds',
        'Input/References/Hrvatin_afteint_refVar.rds'
