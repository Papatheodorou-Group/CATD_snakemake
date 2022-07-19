configfile: 'config.yaml'
include: "Modules/Convert_split/Snakefile"
include: "Modules/Psuedobulk/Snakefile"


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

rule all:
    input:
        outFile
