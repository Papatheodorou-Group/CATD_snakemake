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
1. The cell types should be annotated in the .h5ad or seurat object in the metadata **strictly** under the column name **cellType** (written in camelCase).
2. In the metadata, there needs to be **cellID** and **sampleID** columns, **strictly** under those names. If they don't exist, simply assign rownames to those columns.

### Inputs:
	Input/{sampleName}.h5ad				

	OR				

	Input/{sampleName}.rds    (Should be a Seurat object)

### Outputs:
- Evaluation of selected methods based on selected metrics in config, found in: **Metrics/**
- Plots describing the metrics, found in: **Plots/**
- Resource usage of each step, summarized in: **Benchmarks_summarized.png**
- Individual benchmarks for steps, found in: **Benchmarks/**
- Predictions produced by methods, found in: **Results/**

### Directions

 1. Download the pipeline to your directory of choice.
	`git clone https://github.com/Functional-Genomics/CATD_snakemake.git`

 2. Set up Conda environment with snakemake, use [mamba](https://github.com/mamba-org/mamba) for much faster environment setup.
	  `mamba create -n snakemake snakemake`
 3. Run `basicSetup.sh` to configure conda profile.
 4.  Place the input file into the newly created `Input` directory.
 5. Adjust settings in `config.yaml`.
 6. **(Optional)** Run `getDag.sh` to generate the updated DAG after adjusting config.
 7. Set up cluster resources in `runPip.sh`.
 8. Run the pipeline using `bsub < runPip.sh`.

## Cross-reference
### Description
Uses **two** single-cell references to generate the pseudobulks and references for deconvolution benchmarking.  Important assumptions are
-	The cell types should be annotated in **seurat objects** in the metadata **strictly** under the column name **cellType** (written in camelCase).
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

# Current flow graph
<img src="https://github.com/Functional-Genomics/CATD_snakemake/blob/main/dag.png" alt="drawing">
