# demo_snakemake
A set of scripts to illustrate the use of snakemake


## Purpose

The objectives are:

i) Introducing `snakemake` using the first support:

- https://slides.com/johanneskoester/deck-1
- https://www.slideshare.net/slideshow/introduction-to-snakemake/83052693
- https://moodle.france-bioinformatique.fr/pluginfile.php/346/course/section/47/13_tutoriel_snakemake.html#/

ii) Introducing the CIMENT/GRICAD infrastructure

iii) Put the concepts we've learned into practice through a few use cases.


## Create an account to access to the cluster

ciment_infrastructure.pdf (Figure)

Open a PERSEUS account by clicking on the following link:  

  - https://perseus.univ-grenoble-alpes.fr/create-account/form
  - Choose a login following the recommendations.
  - Use an institutional email address.
  - Choose “formations” as your laboratory.
  - Choose your contract and select a suitable end date.
  - Choose a password.
 
You will receive two emails: 

  - one to join a project, choose the appropriate project,
  - one to validate your email, follow the procedure.
 

## How to start? (180 seconds tutorial)


### Connection to the cluster and job submission

Once your account has become active, log in to the cluster frontend as follows:

```
username=fchdemo
ssh -o ProxyCommand="ssh ${username}@access-gricad.univ-grenoble-alpes.fr -W %h:%p" ${username}@dahu.univ-grenoble-alpes.fr
```

Then, submit an interactive (`-I`) job of 4 cores:

```
oarsub --project groupecalcul -l nodes=1/core=2,walltime=00:03:00 -t fat -I
oarsub --project groupecalcul -t inner=26054906 -l nodes=1/core=8,walltime=00:30:00 -I # mercr.
oarsub --project groupecalcul -t inner=26054914 -l nodes=1/core=8,walltime=00:30:00 -I # jeud.
oarsub --project groupecalcul -t inner=26054926 -l nodes=1/core=4,walltime=04:00:00 -I # vendr.
chandler
oarstat -u cougoulg
```


### Conda environment and working directory

Set up your conda environment as follow:

```
source /home/chuffarf/conda_config.sh
# conda create -n demosnakemake_env
conda activate demosnakemake_env
# mamba install -c anaconda -c bioconda -c conda-forge -c r r-base snakemake=7.32.4 python=3.9 graphviz python-kaleido tenacity plotly 
# pip install smgantt==0.0.5
```

Exercice: Set your own conda environment on your laptop (if needed).

Set up your working directory

```
# Retrieve script and material from github
git clone https://github.com/fchuffar/demo_snakemake.git
# Go to your workdir
cd demo_snakemake
```

### Launch your first workflow

On a worker node or on your laptop:

```
snakemake -s 01st_workflow.py --cores 1 -rp clean
snakemake -s 01st_workflow.py --cores 2 -rpn
snakemake --forceall --dag -s 01st_workflow.py| dot -Tpdf > dag.pdf
smgantt
```

From frontend on many worker nodes:

```
snakemake clean -s 01st_workflow.py --cores 1 -rp
snakemake -s 01st_workflow.py --jobs 50 --cluster "oarsub --project groupecalcul -l nodes=1/core=1,walltime=00:03:00"  --latency-wait 60 -pn
```

Exercice: 

  - Reproduice the *180 seconds tutorial* section.
  - Enhance the *180 seconds tutorial* by adjusting snakemake `cores` argument and `threads` rule parameter, increasing the number of jobs, their duration... 
  - Comment.
  





## Demo (`02nd_worflow.py`)

```
# mamba install -c anaconda -c bioconda -c conda-forge -c r -c brown-data-science r-rmarkdown r-mediation 
snakemake -s 02nd_worflow.py --cores 16 -rpn
snakemake -s 02nd_worflow.py --jobs 50 --cluster "oarsub --project groupecalcul -l nodes=1/core=1,walltime=00:10:00"  --latency-wait 60 -pn
```


### 

