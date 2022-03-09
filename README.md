# EdUSeqAnalysis
 
## clone repository
```
git clone https://github.com/SansamLab/EdUSeqAnalysis.git
```
## run on hpc
```
ml python pandas numpy
sbatch --wrap="snakemake -j 999 --cluster-config config/cluster_config.yml --cluster 'sbatch -A {cluster.account} -p {cluster.partition} --cpus-per-task {cluster.cpus-per-task}  -t {cluster.time} --mem {cluster.mem} --output {cluster.output}'"
```
