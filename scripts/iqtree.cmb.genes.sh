#!/bin/bash
#SBATCH --job-name=iqtree2.glud1
#SBATCH --mem=80G   # Memory per core, use --mem= for memory per node
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20
#SBATCH --time=1-04:00:00   # Use the form DD-HH:MM:SS
#SBATCH --output="%x_%j.o"
#SBATCH --error="%x_%j.e"
#SBATCH --mail-type=ALL

source /homes/lianggao/miniconda3/bin/activate
conda activate iqtree

cd ~/HMW_glutenin/data/phylogenetic_trees/cmb.glud1x.glud1y.partition/
iqtree -s 210604.cmb.glud1x.glud1y.fa -p parti.v1.txt  -m TEST -T 10 -bb 1000 -alrt 1000 > with.partition.out 2> with.partition.err

