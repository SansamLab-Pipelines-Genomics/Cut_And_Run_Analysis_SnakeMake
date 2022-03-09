1) Load slurm and miniconda

```
ml slurm
ml miniconda
```

2A) FIRST TIME ONLY:  Setup environment

```
conda env create -f workflow/envs/qc_trim_align.yml
```

2B) Start environment:

```
conda activate qc_trim_align
```

3)  Modify the config/config.yml

4)  Modify the samples.csv

5) Test with dry run
```
snakemake -npr
```

6) Run with slurm

```
sbatch --wrap="\
snakemake \
-R \
-j 999 \
--cluster-config config/cluster_config.yml \
--cluster '\
sbatch \
-A {cluster.account} \
-p {cluster.partition} \
--cpus-per-task {cluster.cpus-per-task}  \
-t {cluster.time} \
--mem {cluster.mem} \
--output {cluster.output}'"
```

7) When finished, exit environment.

```
conda deactivate
```
