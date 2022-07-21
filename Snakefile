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


outFile = list()
inP = config['Input']['pbulk']
inC = config['Input']['cgen']

if config['seededRun'] == False:
    pbulks = inP + '_pbulks.rds'
    props = inP + '_props.rds'
    c0 = inC + '_C0.rds'
    c1 = inC + '_C1.rds'
    c2 = inC + '_C2.rds'
    refvar = inC + '_refVar.rds'


else:
    pbulks = inP + '_pbulks_seeded.rds'
    props = inP + '_props_seeded.rds'
    c0 = inC + '_C0_seeded.rds'
    c1 = inC + '_C1_seeded.rds'
    c2 = inC + '_C2_seeded.rds'
    refvar = inC + '_refVar_seeded.rds'

outFile.append(pbulks)
outFile.append(props)
outFile.append(c0)
outFile.append(c1)
outFile.append(c2)
outFile.append(refvar)

rule prepModule:
    input:
        outFile,
        "Input/Normalized_tables/Hrvatin_afteint_pbulks_transformed_scaled.rds",
        "Input/Normalized_tables/Hrvatin_afteint_C0_transformed_scaled.rds"

#    output:
#        "passPrep"
#
#    shell:
#        " printf -- 'Delete this file to reactivate the deconvolution preparation module.' > passPrep "
