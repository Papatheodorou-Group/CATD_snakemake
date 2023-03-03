# CATD_snakemake
Critical Assesment of Transcriptomic Deconvonvolution-snakemake pipeline
<pre>

                                                    \`*-.                   
                                                     )  _`-.    Hi          
                                                    .  : `. .               
                                  Hey               : _   '  \              
                     .-=-.          .--.            ; *` _.   `*-._         
         __        .'     '.       /  " )           `-.-'          `-.      
 _     .'  '.     /   .-.   \     /  .-'\             ;       `       `.    
( \   / .-.  \   /   /   \   \   /  /    ^            :.       .        \   
 \ `-` /   \  `-'   /     \   `-`  /                  . \  .   :   .-'   .  
  `-.-`     '.____.'       `.____.'                   '  `+.;  ;  '      :  
                                                      :  '  |    ;       ;-.
https://ascii.co.uk/art/snake                         ; '   : :`-:     _.`* ;
                                                   .*' /  .*' ; .*`- +'  `*'
                                                    `*-*   `*-*  `*-*'     

                                                   https://ascii.co.uk/art/cat
</pre>


Author & maintainer contact : Anna Vathrakokoili Pournara annavp@ebi.ac.uk\
Snakemake Pipeline implementation : Ozgur Beker ozgurbeker@sabanciuniv.edu

If you are interested in our pipeline or you would like to include your method in the CATD pipeline please contact us or create a pull request.

Please use the [issue tracker](https://github.com/Functional-Genomics/CATD_snakemake/issues) if you encounter a problem

# Table of contents
1. [Overview](#overview)
2. [Repo contents](#repocontents)
3. [System requirements](#systemrequirements)
4. [Installation guide](#installationguide)\
 	1.[Instructions](#instructions)\
	2.[Time](#time)
5. [Running the Pipeline](#runningthepipeline)\
 	1.[Self-reference(Task type 1)](#self-reference)\
	2.[Cross-reference(Task type 2)](#cross-reference)\
	3.[Real bulk RNA-seq deconvolution](#realbulk)

## Overview <a name="overview"></a>

The CATD pipeline is a benchmarking pipeline meant to facilitate the assessment of cell-type deconvolution methods(29 now) across different simulation scenarios in a standardised way. It also allows the deconvolution of real bulk samples with many user-unput parameters so that users can deconvolute their own in-house data following our proporsed guidelines.
The pipeline includes:
* Pseudobulk generation methods that allow to create diverse bulk samples and compare deconvolution methods across different scenarios.
* 17 normalization methods implemented in the pipeline for the normalization of the input single-cell reference and the (pseudo)bulk samples
* 4 transformation methods implemented
* 9 DE tests for the selection of marker genes from single-cell reference data(Seurat)
* 29 deconvolution methods
* 7 metrics to assess the results when we test deconvolution methods on pseudo-bulks or when ground truth proportions from real data are available.

For more details check our preprint: 

Vathrakokoili Pournara, A., Miao, Z., Beker, O. Y., Brazma, A. & Papatheodorou, I. Power analysis of cell-type deconvolution methods across tissues. http://biorxiv.org/lookup/doi/10.1101/2023.01.19.523443 (2023) doi:10.1101/2023.01.19.523443.

## Repo contents <a name="repocontents"></a>

The github repository includes: 
* The CATD pipeline with the driver snakefile and main pipeline scripts 

```bash
/bin
.
├── basicSetup.sh  #script for the basic set up 
├── config.yaml #config file with all the user-defined parameters 
├── getDag.sh  #get DAG graphic aka a directed acyclic graph (DAG) of jobs where the edges represent dependencies
├── runPip.sh  #script to run the pipeline
└── Snakefile  #Driver script for the deconvolution assessment pipeline
```
* All the modules used by the pipeline.

Each Module consists of an enviroment file, a snakefile with the rule(s) used in this module & the script to run the method
e.g 


```bash
/EpiDISH
.
├── env.yaml
├── EpiDISH_run.R
└── Snakefile
```

* A README file 

* Demo for running self-reference, cross-reference & real bulk scenarios(tasks)

* Example datasets to run the pipeline


## System Requirements <a name="systemrequirements"></a>

### Hardware requirements
Except for file format conversion and results generation, it is **NOT** recommended to run the pipeline on personal computers. 

**Memory:** File format conversion and results generation parts of the pipeline require relatively less memory (enough to load the dataset into the session, ~4-5 GBs for single cell data of ~30k cells). The memory requirement for pseudo-bulk generation scales with cores provided, so `numThreads*memData` amount of memory should be allocated, with again a good baseline being ~4-5 GBs per thread allocated. For deconvolution with CDSeq, Bisque, CPM, DWLS, MuSiC, SCDC, TIMER, bseqsc and MOMF the memory requirement should be around the same as pseudo-bulk generation. For all other intermediate operations, minimal memory is needed. Note that the memory requirements will definitely change with the size of data, so users are recommended to check logs for specific runs/modules and scale the baseline values given in `config.yaml` accordingly.   

**Threads:** It is highly recommended to run the pipeline with at least 16 threads to avoid long runtimes, especially for pseudo-bulk generation. For bseqsc and CIBERSORT, 3 should be the maximum amount of threads. For DWLS, threads should be equal (ideally) or smaller than the number of cell types. Allocating more threads will not cause any speed-ups for these methods.

For a general comparison of methods and their memory requirements, users can consult Supplementary Figure 9 which was generated with the default values in `config.yaml`. 


### Software requirements

### OS: 
There are no specifications on OS requirements, any OS which can run shell scripts within the pipeline and adhere to remaining requirements should work. The pipeline was tested on Linux with Ubuntu 20.04. 

### Other software:

* git
* conda
* mamba

[Mamba](https://github.com/conda-forge/miniforge#mambaforge) package manager (Mambaforge is recommended) should be installed, with instructions available at the respective repository.

* Snakemake must be installed to run the pipeline, which can be achieved through Mamba with:

	mamba create -n snakemake snakemake

Excluding the initial Snakemake install, conda environments for each module of the pipeline will be handled automatically by Snakemake / Mamba. To achieve this, the setup script will prompt for the absolute path of your `conda.sh` file. The only exception is the `CIBERSORT.R` file which is required for CIBERSORT and bseqsc deconvolution. Users should register and download the code manually, then place it under `Modules/CIBERSORT` and `Modules/bseqsc` respectively.

    
## Installation guide <a name="installationguide"></a>

### Instructions 

 1. Download the pipeline to your directory of choice.

		git clone https://github.com/Functional-Genomics/CATD_snakemake.git

 3. Set up Conda environment with snakemake, pipeline strictly uses [mamba](https://github.com/mamba-org/mamba) for much faster environment setup.

		mamba create -n snakemake snakemake

 4. Build the pipeline.

		make


If you are running the pipeline for the first time you should use one sample to run in order for all the environments to be installed fast( each module/method in our pipeline has its own environment. 

 5.  Place the input file in the `Input` directory that has been created ( example input file under example/first_input.h5ad  )
 6. Adjust settings in `config.yaml`. For the first run use the default settings that we have there already
 7. **(Optional)** Run `getDag.sh` to generate the updated DAG after adjusting config. 
 8. **(Optional)** If on cluster set up cluster profile if you haven't, instructions available [here](https://github.com/Snakemake-Profiles/lsf).
 9. Run the pipeline using `bsub < runPip.sh` or through `snakemake --use-conda --cores [N]` if on local.


### Time
**Installation time & first mock run to install environments** :
The steps(1-2-3-4) should take few seconds to complete./
If all the environments are included in the first run it will take around 25-30 minutes, roughly 1min/environment is a good estimate.(step 9)


# Running the pipeline <a name="runningthepipeline"></a>

**IMPORTANT**: As mentioned before if running for the first time, use one sample only as the environments are
installed throughout the workflow.(See Installation guide)


## Self-reference deconvolution
### Description
Uses **one** single-cell reference to generate both the pseudobulks and references and aims to aid the benchmark of deconvolution methods./
Important assumptions:
*  The cell types should be annotated in the .h5ad or seurat object in the metadata **strictly** under the column name **cellType** (written in camelCase).
*  In the metadata, there needs to be **cellID** and **sampleID** columns, **strictly** under those names. If they don't exist, simply assign rownames to those columns.

### Inputs:
	Input/{sampleName}.h5ad				

	OR				

	Input/{sampleName}_seurat.rds

### Outputs:
- Evaluation of selected methods based on selected metrics in config, found in: **Metrics/**
- Plots describing the metrics, found in: **Plots/**
- Resource usage of each step, summarized in: **{sampleName}_benchmarks_summarized.png**
- Individual benchmarks for steps, found in: **Benchmarks/**
- Predictions produced by methods, found in: **Results/**


## Cross-reference deconvolution

### Description
Uses **two** single-cell references to generate the pseudobulks and references for deconvolution benchmarking. Aims to study the batch effects caused because of differences between the bulk and the reference(technology effect,sample and study effect)./

Important assumptions are:/
-	All assumptions in the self-reference part/
-	The **levels** (i.e unique list) of cell types must be the **same** in both references provided.

### Inputs:
	Input/Cell_splits/{sampleName}_gen.rds		(Will be used to generate psuedobulks)				
	Input/Cell_splits/{sampleName}_C0.rds    	(Will be used to generate references)

### Outputs:
Same as self-reference.

### Instructions
Same as self-reference, except after the **3rd** step, note the following directory:

	Input/Cell_splits

**The input files should go in this folder**. Make sure that the inputs **conform to the standards written in the 'Inputs' section above**. Then continue with the **5th** step.

## Real Bulk

### Description
Uses **one** reference single cell matrix with **user-defined** bulks and (optional) **known** proportions for deconvolution benchmarking. 
Assumptions are:/
- All assumptions in the self reference part/
- The **rownames** (cell types) in the proportions should  be the **same** as the cell types annotated in the reference
### Inputs
	Input/{sampleName}.h5ad  /   Input/{sampleName}_seurat.rds 

	OR

	Input/Cell_splits/{sampleName}_C0.rds

	AND

	Input/Psuedobulks/{sampleName}_pbulks.rds
	Input/Psuedobulks/{sampleName}_props.rds


The first two options will only use **half** of the data to generate references. The third will use **all** of the data to generate the reference. Alongside the reference, you need to input the pseudo-bulks inside the folder specified under those names.

### Outputs
Same as self-reference.

### Instructions
Same as self-reference, except after the **3rd** step, note the following directory:

	Input/Psuedobulks

This is where the bulks should go. Also, **enable realBulk in config.** If you wish to use all the data for the reference, go to:

	Input/Cell_splits

and place the reference here under the name **{sampleName}_C0.rds**. Otherwise, you can use the `Input` directory with **{sampleName}.h5ad / {sampleName}_seurat.rds**


# Appendices
## Running CIBERSORT and bseqsc in the pipeline
Users need to download the CIBERSORT.R code manually and place it within the modules (`Modules/CIBERSORT` and `Modules/bseqsc`) to run these methods.

## Sample directed acyclic graph (DAG) of the pipeline
The image below is generated by selecting only a couple of methods. If all methods are selected, snakemake cannot generate a pretty DAG as all deconvolution methods run in parallel.
<img src="https://github.com/Functional-Genomics/CATD_snakemake/blob/main/dag.png"/>

## Comprehensive Joblist
The full list of jobs looks like this:

	job                           count    min threads    max threads
	--------------------------  -------  -------------  -------------
	Autogenes_run                     1              1              1
	Bisque_run                        1              1              1
	CDSeq_run                         1             16             16
	CIBERSORT_run                     1              3              3
	CPM_run                           1             16             16
	CellmixInstall                    1              1              1
	DCQ_run                           1              1              1
	DSA_run                           1              1              1
	DWLS_run                          1             16             16
	DeconRNASeq_run                   1              1              1
	EPIC_run                          1              1              1
	EpiDISH_run                       1              1              1
	FARDEEP_run                       1              1              1
	Momf_run                          1              1              1
	MuSiC_run                         1              1              1
	NNLS_run                          1              1              1
	OLS_run                           1              1              1
	OmniInstall                       1              1              1
	RLR_run                           1              1              1
	SCDC_run                          1              1              1
	TIMER_run                         1              1              1
	all                               1              1              1
	bseqsc_run                        1              3              3
	convertAnndata                    1              1              1
	debCAM_C1                         1              1              1
	debCAM_marker                     1              1              1
	deconf_run                        1              1              1
	dtangle_run                       1              1              1
	elasticNET_run                    1              1              1
	exploreBenchmarks                 1              1              1
	exploreResults                    1              1              1
	generatePsuedobulks               1             32             32
	generateReferences                1             32             32
	lasso_run                         1              1              1
	proportionsInAdmixture_run        1              1              1
	ridge_run                         1              1              1
	scaleTransform_T_C                1             32             32
	splitCells                        1              1              1
	ssFrobenius_run                   1              1              1
	ssKL_run                          1              1              1
	sumRun                            1              1              1
	vioplotResults                    1              1              1
	visualizeResults                  1              1              1
	total                            43              1             32

## Cleaning all outputs
Use `make clean` to delete all outputs from a workflow. Keep in mind that this will delete **ALL** outputs including results.
