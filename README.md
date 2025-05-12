## Dog genome-wide association study and heritability workflow (dog-gwas-heritability-nextflow) 
Nextflow workflow to conduct genome-wide association studies in dogs using a mixed linear model with leave-one-chromosome-out approach (GCTA MLMA LOCO) and heritability analyses using restricted maximum likelihood

### Input files 

The input files include the following:

  * A phenotype file containing the column `dog_id` and columns for each phenotype of interest to be analyzed (tab-delimited format with headers: dog_id "\t" phenoA "\t" phenoB)

```
dog_id  F1      F2      F3      F4      F5      F6      F7      F8
28      5.114824997780695       4.11126844119854        -1.367657516345442      3.2407159840647464      15.110542318327791      3.252383263718634       -8.958823306743273      -0.8692999350381554
46      12.041305639988916      10.18177012301823       -1.415599633101153      0.5474706316540905      17.082930711862293      9.72219562232409        -11.800388098421916     -5.28880991642094
95      9.786001103329315       16.189784949097383      -1.910598665812524      -11.137085834669682     23.089078417388883      -1.8888784705935753     -10.38988615567527      -5.28880991642094
98      1.884329900709555       NA      NA      NA      NA      NA      NA      NA
100     11.911674129898518      10.67338676061298       -1.910598665812524      -7.498135673754153      17.65724309943409       6.671555009460411       -13.032302510030297     -2.381005844829711


```
  * A quantitative covariate file containing the columns `IID`, `FID`, and columns for each quantitative covariate of interest (such as age and weight) but without header in tab-delimited format
    
```
0       28      1.68407 3
0       46      1.5219800000000001      3
0       95      7.53571 3000000001      3
0       98      13.8187 2
```
  * A discrete covariate file containing the columns `IID`, `FID`, and columns for each discrete covariate of interest (such as sex) but without header in tab-delimited format

```
0       28      female  yes
0       46      male    yes
0       95      female  yes
0       98      female  yes

```
  * Genetic plinkset (filtering of plinkset occurs within the workflow)
  * Annotation file in the following format (tab-delimited file containing chr, start position, end position, gene name)

 ```
1       13565063        13577444        SERPINB13
1       13587694        13611747        SERPINB12
1       13639168        13716725        SERPINB5
1       13671047        13674721        LOC119870705
1       13716837        13751163        VPS4B
1       13766594        13804675        KDSR
1       13811586        13978482        BCL2


```

### Parameters 
Default values are listed that can be changed upon submission of workflow, and path to input files are required for --geneticset, --annotation, --pheno_dateset, --covar, and --qcovar

```
params.geneticset = "/scratch/vsohrab/darwins_dogs/darwins_dogs_genetic_set_2024/DarwinsDogs_2024_N-3277_canfam4_gp-0.70_biallelic.{bed,bim,fam}"

params.maf = 0.01 (minor allele frequency threshold - default value is 0.01)

params.geno = 0.05 (genotyping threshold - default value is 0.05)

params.hwe = 1e-20 (Hardy-Weinberg equilibrium p-value - default value is 1e-20)

params.autosome_num = 38 (number of autosomes - default value for dogs is 38)

params.ld_window_kb = 250 (specifying the length of the region for segment-based LD score calucation of GCTA - default value is 250 kb) 

params.clump_pval = 1e-6 (significance threshold for index SNPs - default value is 1e-6)

params.clump_kb = 250 (physical distance threshold for clumping
 - default value is 250 kb)

params.clump_rsquared = 0.2 (LD threshold for clumping - default is 0.2)

params.annotation = "/scratch/vsohrab/reference/UU_Cfam_GSD_1.0_ROSY.refSeq.ensformat.genes.validchr.bed"

params.pheno_dataset = /scratch/vsohrab/dog_gwas/DarwinsArk_8Factors_N-3277.tsv

params.covar = /scratch/vsohrab/dog_gwas/DarwinsArk_AgeComposed_Height_N-3277.tsv

params.qcovar = /scratch/vsohrab/dog_gwas/DarwinsArk_Sex_NeuterStatus_N-3277.tsv 

```

#### Note
Each workflow run is unique to a specific discrete covariate and quantitative covariate set pair (for example, if I have 10 phenotypes to analyze where quantitative covariates are age and weight and discrete covariate is sex; then, I can include those 10 phenotypes in a single phenotype input file with dog_id, pheno1, pheno2, pheno3,...,pheno10; however, if I have another set of phenotypes with a different set of covariates, then I would need to run those separately. If quantitative covariate is age and discrete covariate is sex, then I will create a separate phenotype file for that workflow run with dog_id, phenoA, phenoB). If running the same phenotype with and/or without a particular covariate, embed the covar and qcovar name in the phenotype input file column name, because the column names in the input phenotype file is used to create the output files (and if those are the same names, then the will be overwritten with the most recent workflow run when running the same phenotype with different set of covariates in the same working directory)


### Output files

* GWAS output file(s) are reported in *association* folder
* Heritability output file(s) are reported in *reml* folder
* Phenotype file(s) generated by workflow can be found in *phenofiles* folder
* Clump output file(s) are reported in *clump_pval-$pval_kb-$kb_r2-$r2* folder. For example if all default parameters are used, clump output files can be found in *clump_pval-0.000001_kb-250_r2-0.2*
* Manhattan plots and qqplots are generated in *plots* folder
* For viewing an example log file from GWAS result, please check the *logs* folder (ie to view sample size of GWAS analysis for a particular phenotype)


### Example code to run genetic analysis: 

```
/PATH/TO/nextflow run /PATH/TO/dog_gwas_heritability_multiple_phenotypes.nf \
-c /PATH/TO/nextflow_asu_dog_gwas_heritability_metabolites.config \
-resume -w /PATH/TO/work \
--geneticset "/PATH/TO/genetic_set_prefix.{bed,bim,fam}" \
--pheno_dataset /PATH/TO/phenotype_input_file.tsv \
--covar /PATH/TO/discrete_covar_input_file.tsv \
--qcovar /PATH/TO/quant_covar_input_file.tsv \
--annotation /PATH/TO/annotation.bed \
--outdir /PATH/TO/workdir \
--maf ${maf} --geno ${geno} --hwe ${hwe}
```

