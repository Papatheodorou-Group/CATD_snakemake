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

Maintainer contact : Anna Vathrakokoili Pournara annavp@ebi.ac.uk\
Software Authors ✍🏼 ✍🏼: Chichao Miao, Anna Vathrakokoili Pournara & Ozgur Beker\
Snakemake Pipeline implementation 👨‍💻 : Ozgur Beker ozgurbeker@sabanciuniv.edu

Many Thanks to Nadja Nolte for testing and further developping.


If you are interested in our pipeline or you would like to include your method in the CATD pipeline please contact us or create a pull request.

Please use the [issue tracker](https://github.com/Functional-Genomics/CATD_snakemake/issues) if you encounter a problem

# Table of contents
 * [Overview](#overview)
 * [Repo contents](#repocontents)
 * [System requirements](#systemrequirements)
 * [Installation guide](#installationguide)
 	* [Instructions](#instructions)
	* [Time](#time)
 * [Running the pipeline](#runningthepipeline)
 	* [Self-reference deconvolution](#self-referencedeconvolution)
	* [Cross-reference deconvolution](#cross-referencedeconvolution)
	* [Real bulk RNA-seq deconvolution](#realbulkRNA-seqdeconvolution)
 * [Appendices](#appendices)

## Overview <a name="overview"></a>

The CATD pipeline is a benchmarking pipeline meant to facilitate the assessment of cell-type deconvolution methods(29 now) across different simulation scenarios in a standardised way. It also allows the deconvolution of real bulk samples with many user-unput parameters so that users can deconvolute their own in-house data following our proporsed guidelines.
The pipeline includes:
* Pseudobulk generation methods that allow to create diverse bulk samples and compare deconvolution methods across different scenarios.
* 17 normalization methods implemented in the pipeline for the normalization of the input single-cell reference and the (pseudo)bulk samples
* 4 transformation methods implemented
* 9 DE tests for the selection of marker genes from single-cell reference data(Seurat)
* 29 deconvolution methods
* 7 metrics to assess the results when we test deconvolution methods on pseudo-bulks or when ground truth proportions from real data are available.

For more details check our [preprint](https://www.biorxiv.org/content/10.1101/2023.01.19.523443v1): 

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
/Modules/EpiDISH
.
├── env.yaml
├── EpiDISH_run.R
└── Snakefile
```

* A README file containing between other things : 

	* Instructions/Demo for running self-reference, cross-reference & real bulk scenarios(tasks)

	* examples datasets deposited on google drive


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

[Mamba](https://github.com/conda-forge/miniforge#mambaforge) package manager (Mambaforge is recommended) should be installed, with instructions available at the respective repository.

* Snakemake must be installed to run the pipeline, which can be achieved through Mamba with:

	mamba create -n snakemake snakemake

Excluding the initial Snakemake install, conda environments for each module of the pipeline will be handled automatically by Snakemake / Mamba. To achieve this, the setup script will prompt for the absolute path of your `conda.sh` file. The only exception is the `CIBERSORT.R` file which is required for CIBERSORT and bseqsc deconvolution. Users should register and download the code manually, then place it under `Modules/CIBERSORT` and `Modules/bseqsc` respectively.

    
## Installation guide <a name="installationguide"></a>

### Instructions 

 1. Download the pipeline to your directory of choice.

		git clone https://github.com/Functional-Genomics/CATD_snakemake.git

 2. Set up Conda environment with snakemake, pipeline strictly uses [mamba](https://github.com/mamba-org/mamba) for much faster environment setup.

		mamba create -n snakemake snakemake

 3. Build the pipeline.

		make


If you are running the pipeline for the first time you should use one sample to run in order for all the environments to be installed fast( each module/method in our pipeline has its own environment. 

 4.  Place the input file in the `Input` directory that has been created (follow the [link](https://drive.google.com/drive/u/0/folders/14rCp3kxUsT3J4DIDg4mW-1QhlcAikE88) to the example file - use either XinY_4Decon.h5ad or XinY_4Decon_seurat.rds )
 5. Adjust settings in `config.yaml`like suggested: (familiarise yourself with the config.yaml file and all the parameters that you can choose)
 Specifically:
 	**sampleNames**: #Enter the part BEFORE file extension (ex. 'data.h5ad' should be entered as 'data')
		- XinY_4Decon # in this case
        **pbulkParam:**
  pbulkMode: 2  #Different mode for building pbulks
  cellCount: 100 #How many cells should be used to create a bulk sample
  nSamples: 1 #How many different samples should be generated #first time only generate one sample
  propVar: 500 #Variances inbetween proportions, check vignette (EFFECTIVE IN MODE 2 ONLY, enter negative value to switch to min max prop mode)
  sampleCT: 0 #Sampling for cell types allowed or not (1 or 0) (EFFECTIVE IN MODE 2 ONLY)
	
 6. **(Optional)** Run `getDag.sh` to generate the updated DAG after adjusting config. 
 	./getDag.sh 
	this will output a dag.png where you can observe the jobs that will be executed
 7. **(Optional)** If on cluster set up cluster profile if you haven't, instructions available [here](https://github.com/Snakemake-Profiles/lsf).
 8. **(Optional)** Dry run: 
 	conda activate snakemake 
	snakemake -n
 9. Run the pipeline using `bsub < runPip.sh` (LSF) or through `snakemake --use-conda --cores [N]` if on local.


### Time <a name="Time"></a>
The steps(1-2-3) should take few seconds to complete.
If all the environments are included in the first run it will take around 25-30 minutes,roughly 1min/environment is a good estimate.(step 9)
When the job(step 9) finishes all the environments will have been installed and users are ready to run actual experiments and test the pipeline

# Running the pipeline <a name="runningthepipeline"></a>
**IMPORTANT**:As mentioned above, if running the pipeline for the first time, use one sample only as all the environments are
installed throughout the workflow.


## Self-reference deconvolution <a name="self-referencedeconvolution"></a>

### Description
Uses **one** single-cell reference to generate both the pseudobulks and the reference for deconvolution benchmarking. Important assumptions:
*  The cell types should be annotated in the .h5ad or seurat object in the metadata **strictly** under the column name **cellType** (written in camelCase).
*  In the metadata, there needs to be **cellID** and **sampleID** columns, **strictly** under those names. If they don't exist, simply assign rownames to those columns.

### Inputs:
	Input/{sampleName}.h5ad				

	OR				

	Input/{sampleName}_seurat.rds

###Example data:

[XinY_4Decon.h5ad]()
[XinY_4Decon_seurat.rds]()

### Outputs:
- Evaluation of selected methods based on selected metrics in config, found in: **Metrics/**
- Plots describing the metrics, found in: **Plots/**
- Resource usage of each step, summarized in: **{sampleName}_benchmarks_summarized.png**
- Individual benchmarks for steps, found in: **Benchmarks/**
- Predictions produced by methods, found in: **Results/**

### Instructions

-Start from step **4** from above
-Edit again the 'config.yaml' file :
	This time you can change/tune the parameters e.g increase the **nSamples** to 100 or more , increase **cellCount** to 1000 etc
	based on the pseudobulks you want to create, you can also edit tranformation/normalization/deconvolution methods parameters
-step **6** and **8** from above are optional again but recommended.
-run the self-reference task: step **9**


## Cross-reference deconvolution <a name="cross-referencedeconvolution"></a>

### Description
Uses **two** single-cell references to generate the pseudobulks and references for deconvolution benchmarking.  Important assumptions are:
-	All assumptions in the self-reference part
-	The **levels** (i.e unique list) of cell types must be the **same** in both references provided.

### Inputs:
	Input/Cell_splits/{sampleName}_gen.rds		(Will be used to generate psuedobulks)				
	Input/Cell_splits/{sampleName}_C0.rds    	(Will be used to generate references)

###Example data:

[BaronHuman_Sege_cross_gen.rds](https://drive.google.com/drive/u/0/folders/1mQAsxvQywW3Qt4rm1_SuNADOk5ylJCjS) #Segerstope dataset from where pseudobulks will be created
[BaronHuman_Sege_cross_C0.rds](https://drive.google.com/drive/u/0/folders/1mQAsxvQywW3Qt4rm1_SuNADOk5ylJCjS) #BaronHuman dataset from where the reference will be generated

### Outputs:
Same as self-reference.

### Instructions
Same as self-reference, except after the **4** step, add the input files in the following directory (instead of just Input/):

	Input/Cell_splits

**The input files should go in the above folder**. Make sure that the inputs **conform to the standards written in the 'Inputs' section above**. Then continue with the **6th** step.





## Real bulk RNA-Seq deconvolution <a name="realbulkrna-seqdeconvolution"></a>

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


# Appendices <a name="appendices"></a>
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
