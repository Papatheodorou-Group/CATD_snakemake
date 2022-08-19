
# CATD_snakemake
<pre>

                                                    \`*-.                   
                                                     )  _`-.                
                                                    .  : `. .               
                                                    : _   '  \              
                     .-=-.          .--.            ; *` _.   `*-._         
         __        .'     '.       /  " )           `-.-'          `-.      
 _     .'  '.     /   .-.   \     /  .-'\             ;       `       `.    
( \   / .-.  \   /   /   \   \   /  /    ^            :.       .        \   
 \ `-` /   \  `-'   /     \   `-`  /                  . \  .   :   .-'   .  
  `-.-`     '.____.'       `.____.'                   '  `+.;  ;  '      :  
                                                      :  '  |    ;       ;-.
                                                      ; '   : :`-:     _.`* ;
                                                   .*' /  .*' ; .*`- +'  `*'
                                                    `*-*   `*-*  `*-*'    



Credits:
snake: jgs
cat: bug


</pre>

# Running the pipeline
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
- Resource usage of each step, summarized in: **Benchmarks_summarized.png**
- Individual benchmarks for steps, found in: **Benchmarks/**
- Predictions produced by methods, found in: **Results/**

### Directions

 1. Download the pipeline to your directory of choice.

	    git clone https://github.com/Functional-Genomics/CATD_snakemake.git

 3. Set up Conda environment with snakemake, pipeline strictly uses [mamba](https://github.com/mamba-org/mamba) for much faster environment setup.

		  mamba create -n snakemake snakemake

 4. Run `basicSetup.sh` to configure conda profile.
 5.  Place the input file into the newly created `Input` directory.
 6. Adjust settings in `config.yaml`.
 7. **(Optional)** Run `getDag.sh` to generate the updated DAG after adjusting config.
 8. Set up cluster resources in `runPip.sh`.
 9. Run the pipeline using `bsub < runPip.sh`.

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

### Directions
Same as self-reference, except after the **3rd** step, create a directory named `Cell_splits` within input, using:

	mkdir -p Input/Cell_splits

Then, **place the input files in this folder**. Make sure that the inputs **conform to the standards written in the 'Inputs' section above**. Then continue with the **5th** step.

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


The first two options will only use **half** of the data to generate references. The third will use **all** of the data to generate the reference. Alongside the reference, you need to input the pseudo-bulks inside the folder specified under those names (note that there is a typo in the "Psuedobulks", still the folder should be under that name.)

### Outputs
Same as self-reference.

### Directions
Same as self-reference, except after the **3rd** step, create a directory named `Psuedobulks` within input, using:

	mkdir -p Input/Psuedobulks

This is where the bulks should go. If you wish to use all the data for the reference, also create a directory named Cell_splits using

	mkdir -p Input/Cell_splits

and place the reference here under the name **{sampleName}_C0.rds**
