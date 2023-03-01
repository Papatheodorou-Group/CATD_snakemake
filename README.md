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
                                                    `*-*   `*-*  `*-*'     https://ascii.co.uk/art/bug

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
5. [Functionalities](#functionalities)\
 	1.[Self-reference(Task type 1)](#self-reference)\
	2.[Cross-reference(Task type 2)](#cross-reference)\
	3.[Real bulk RNA-seq deconvolution](#realbulk)

## Overview <a name="overview"></a>

The CATD pipeline is a benchmarking pipeline meant to facilitate the assessment of cell-type deconvolution methods(29 now) across different simulation scenarios in a standardised way. It also allows the deconvolution of real bulk samples with many user-unput parameters so that users can deconvolute their own in-house data following our proporsed guidelines.
The pipeline includes:
* Pseudobulk generation methods that allow to create diverse bulk samples and compare deconvolution methods across different scenarios.
* Seventeen(17) normalization methods implemented in the pipeline for the normalization of the input single-cell reference and the (pseudo)bulk samples
* Four(4) transformation methods implemented
* Nine(9) DE tests for the selection of marker genes from single-cell reference data(Seurat)
* Twenty-nine(29) deconvolution methods
* Seven(7) metrics to assess the results when we test deconvolution methods on pseudo-bulks or when ground truth proportions from real data are available.

For more details check our preprint: 

Vathrakokoili Pournara, A., Miao, Z., Beker, O. Y., Brazma, A. & Papatheodorou, I. Power analysis of cell-type deconvolution methods across tissues. http://biorxiv.org/lookup/doi/10.1101/2023.01.19.523443 (2023) doi:10.1101/2023.01.19.523443.

## Repo contents <a name="repocontents"></a>

Describe the folders in this repo

## System Requirements <a name="systemrequirements"></a>

Hardware requirements

Software requirements
    
    
## Installation guide <a name="installationguide"></a>

 1. Download the pipeline to your directory of choice.

		git clone https://github.com/Functional-Genomics/CATD_snakemake.git

 3. Set up Conda environment with snakemake, pipeline strictly uses [mamba](https://github.com/mamba-org/mamba) for much faster environment setup.

		mamba create -n snakemake snakemake

 4. Build the pipeline.

		make

 5.  Place the input file in the `Input` directory.
 6. Adjust settings in `config.yaml`.
 7. **(Optional)** Run `getDag.sh` to generate the updated DAG after adjusting config.
 8. **(Optional)** If on cluster set up cluster profile if you haven't, instructions available [here](https://github.com/Snakemake-Profiles/lsf).
 9. Run the pipeline using `bsub < runPip.sh` or through `snakemake --use-conda --cores [N]` if on local.



## Functionalities <a name="functionalities"></a>


# Running the pipeline
**IMPORTANT**: If running for the first time, use one sample only as the environments are
installed throughout the workflow.


## Self-reference
### Description
Uses **one** single-cell reference to generate the pseudobulks and references for deconvolution benchmarking. Important assumptions:
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


## Cross-reference

### Description
Uses **two** single-cell references to generate the pseudobulks and references for deconvolution benchmarking.  Important assumptions are:
-	All assumptions in the self-reference part
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
Uses **one** reference single cell matrix with **user-defined** bulks and **known** proportions for deconvolution benchmarking. Assumptions are:
- All assumptions in the self reference part
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
