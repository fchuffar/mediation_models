# mediation_models

A set of scripts to illustrate mediation models


 

## How to start? 


### Conda environment and working directory

Set up your conda environment as follow:

```
# conda create -n mediationmodels_env
conda activate mediationmodels_env
# mamba install -c anaconda -c bioconda -c conda-forge -c r  -c brown-data-science r-base r-rmarkdown r-mediation snakemake=7.32.4 python=3.9 graphviz python-kaleido tenacity plotly 
# pip install smgantt==0.0.5
```


Set up your working directory

```
# Retrieve script and material from github
git clone https://github.com/fchuffar/mediation_models.git
# Go to your workdir
cd mediation_models
```

### 1st simulation


```
snakemake -s 01st_simu.py --cores 16 -rp clean
snakemake -s 01st_simu.py --cores 16 -rpn
# snakemake -s 01st_simu.py --jobs 50 --cluster "oarsub --project epimed -l nodes=1/core=1,walltime=00:10:00"  --latency-wait 60 -pn
```

Exercice:

Targeted dataset is `datasets/sim_res_2000_1_100_1_500_1_50_1_rep1.txt`

  - Launch `mediation_analysis_tp.Rmd` (`rmarkdown::render("mediation_analysis_tp.Rmd")`)
  - Study the nature of the variables CTCF, c3D and mRNA
  - Adapt the `mediation_analysis_tp.Rmd` script and propose a `mediation_analysis_glm.Rmd` script better able to handle this type of data.
  - What is the impact of the changes on the effects of the various mediation scheme models and on the ADE, as well as on the ACME?


### 2nd simulation

Study the datset `datasets/sim_res_2000_-25_100_0_500_1_50_1_rep0.txt`

```
snakemake -s 02nd_simu.py --cores 16 -rp clean
snakemake -s 02nd_simu.py --cores 16 -rpn
```

### 3rd simulation

```
snakemake -s 03rd_simu.py --cores 16 -rp clean
snakemake -s 03rd_simu.py --cores 16 -rpn
cat sim_res_2000_-25_100_0_2000_1_50_1_rep0.txt > datasets/sim_res_3.txt ; cat sim_res_2000_1_100_1_0_1_50_1_rep0.txt >> datasets/sim_res_3.txt
```
