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
include: "Modules/Scale_transform/Snakefile"
include: "Modules/debCAM/Snakefile"


outFile = list()
inP = "Input/Psuedobulks/" + config['sampleName']
inC = "Input/References/" + config['sampleName']
inN = "Input/Normalized_tables/" + config['sampleName']


if not config['seededRun']:
    props = inP + '_props.rds'
    c1 = inC + '_C1.rds'
    c2 = inC + '_C2.rds'
    refvar = inC + '_refVar.rds'
    c0 = inN + '_C0'
    pbulks = inN + '_pbulks'

    if config["stParam"]["scaleFirst"]:
        c0 = c0 + '_scaled_transformed.rds'
        pbulks = pbulks + '_scaled_transformed.rds'

    else:
        c0 = c0 + '_transformed_scaled.rds'
        pbulks = pbulks + '_transformed_scaled.rds'


else:
    props = inP + '_props_seeded.rds'
    c1 = inC + '_C1_seeded.rds'
    c2 = inC + '_C2_seeded.rds'
    refvar = inC + '_refVar_seeded.rds'
    c0 = inN + '_C0_seeded'
    pbulks = inN + '_pbulks_seeded'

    if config["stParam"]["scaleFirst"]:
        c0 = c0 + '_scaled_transformed.rds'
        pbulks = pbulks + '_scaled_transformed.rds'

    else:
        c0 = c0 + '_transformed_scaled.rds'
        pbulks = pbulks + '_transformed_scaled.rds'


outFile.append(pbulks)
outFile.append(props)
outFile.append(c0)
outFile.append(c1)
outFile.append(c2)
outFile.append(refvar)


rule all:
    input:
        #outFile,
        "Output/Hrvatin_afteint_debCAM_unsupervised.txt"

#    output:
#        "passPrep"
#
#    shell:
#        " printf -- 'Delete this file to reactivate the deconvolution preparation module.' > passPrep "
