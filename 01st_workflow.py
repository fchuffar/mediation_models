rule target:
  threads: 1
  input:
    "sim_results_01.txt",
    "sim_results_02.txt",
    "sim_results_03.txt",
    # [f"sim_results_{id}.txt" for id in range(40)]
  shell:"""
echo "Tout est accompli." 
"""

rule clean:
  threads: 1
  shell:"rm -Rf *.txt"

rule run_simulation:
  output:"sim_results_{id}.txt",
  threads: 1
  shell:"""
sleep 3
echo $HOSTNAME > {output}
"""

localrules: target clean
