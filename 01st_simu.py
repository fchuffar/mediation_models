localrules: target, clean

import os 
import os.path
prefix = os.getcwd()

# parameter sweeps
ncells = ["100000"]
ncells = ["500"]
# ncells = ["100"]
mis    = ["1"]
sijs   = ["100"] 
ds     = ["1"]
fs     = ["50"]
lambs  = ["1"]
gamma1s = ["0", "1"]
gamma2s = ["0", "500"]
seeds   = ["0", "1", "2"]
# seeds = [f"{seed}" for seed in range(8)]
htmls = [f"{prefix}/mediation_analysis_lm_{ncell}_{i}_{sij}_{gamma1}_{gamma2}_{d}_{f}_{lamb}_rep{seed}.html" for ncell in ncells for i in mis for sij in sijs for gamma1 in gamma1s for gamma2 in gamma2s for d in ds for f in fs for lamb in lambs for seed in seeds]

# Go snakemake, GO!!
rule target:
  threads: 1
  message: "-- Execute the DAG. --"
  input:
    htmls
  shell:"""
export OMP_NUM_THREADS=1
source /home/chuffarf/conda_config.sh
conda activate demosnakemake_env

RCODE="rmarkdown::render('metaanalysis.Rmd')"
echo $RCODE | Rscript -
snakemake --forceall --dag -s 01st_simu.py | dot -Tpdf > dag.pdf
smgantt
echo "done." 
"""

rule run_simulation:
    input:
      sim   = "{prefix}/simulator.py",
    output:
      sres  = "{prefix}/sim_res_{ncell}_{i}_{sij}_{gamma1}_{gamma2}_{d}_{f}_{lamb}_rep{seed}.txt",
    log:
      "{prefix}/sim_log_{ncell}_{i}_{sij}_{gamma1}_{gamma2}_{d}_{f}_{lamb}_rep{seed}.log"
    threads: 1
    shell:"""
export OMP_NUM_THREADS=1
source /home/chuffarf/conda_config.sh
conda activate demosnakemake_env

cd {wildcards.prefix}

RCODE="
  ncell = '{wildcards.ncell}'; i = '{wildcards.i}'; sij = '{wildcards.sij}'; gamma1 = '{wildcards.gamma1}'; gamma2 = '{wildcards.gamma2}'; d = '{wildcards.d}'; f = '{wildcards.f}'; lamb = '{wildcards.lamb}'; seed = {wildcards.seed}; source('params.R');
  rmarkdown::render('run_simulation.Rmd', output_file=paste0('run_simulation_', suffix, '.html'), intermediates_dir=paste0('tmp_', suffix))"
echo $RCODE | Rscript - 2>&1 > {log}
"""

rule mediation_analysis_lm:
    input:
      rmd   = "{prefix}/mediation_analysis_lm.Rmd",
      par   = "{prefix}/params.R",
      sres  = "{prefix}/sim_res_{ncell}_{i}_{sij}_{gamma1}_{gamma2}_{d}_{f}_{lamb}_rep{seed}.txt",
    output:
      html  = "{prefix}/mediation_analysis_lm_{ncell}_{i}_{sij}_{gamma1}_{gamma2}_{d}_{f}_{lamb}_rep{seed}.html"    ,
      rout  = "{prefix}/mediation_analysis_lm_{ncell}_{i}_{sij}_{gamma1}_{gamma2}_{d}_{f}_{lamb}_rep{seed}.Rout"    ,
      infos = "{prefix}/infos_{ncell}_{i}_{sij}_{gamma1}_{gamma2}_{d}_{f}_{lamb}_rep{seed}.rds",
      fit   = "{prefix}/fitmed_{ncell}_{i}_{sij}_{gamma1}_{gamma2}_{d}_{f}_{lamb}_rep{seed}.rds",
    log:
      rout  = "{prefix}/mediation_analysis_lm_{ncell}_{i}_{sij}_{gamma1}_{gamma2}_{d}_{f}_{lamb}_rep{seed}.Rout"    ,
    threads: 1
    shell:"""
export OMP_NUM_THREADS=1
source /home/chuffarf/conda_config.sh
conda activate demosnakemake_env

cd {wildcards.prefix}

RCODE="
  ncell = '{wildcards.ncell}'; i = '{wildcards.i}'; sij = '{wildcards.sij}'; gamma1 = '{wildcards.gamma1}'; gamma2 = '{wildcards.gamma2}'; d = '{wildcards.d}'; f = '{wildcards.f}'; lamb = '{wildcards.lamb}'; seed = {wildcards.seed}; source('params.R');
  rmarkdown::render('mediation_analysis_lm.Rmd', output_file=paste0('mediation_analysis_lm_', suffix, '.html'), intermediates_dir=paste0('tmp_', suffix))"
echo $RCODE | Rscript - 2>&1 > {log}
"""

rule clean:
  threads: 1
  shell:"rm -Rf *.txt *.rds *.Rout *.html *.log tmp_* .snakemake"
