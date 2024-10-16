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

  - Study the nature of the variables X, M and Y
  - Adapt the `mediation_analysis_lm.Rmd` script and propose a `mediation_analysis_glm.Rmd` script better able to handle this type of data.
  - What is the impact of the changes on the effects of the various mediation scheme models and on the ADE, as well as on the ACME?

### 2nd simulation


```
snakemake -s 02nd_simu.py --cores 16 -rp clean
snakemake -s 02nd_simu.py --cores 16 -rpn
```

