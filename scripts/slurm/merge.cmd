#!/bin/bash
#SBATCH -N 1 # node count
#SBATCH --ntasks-per-node=1
#SBATCH -t 239:59:00
#SBATCH --mem=120GB
# sends mail when process begins, and 
# when it ends. Make sure you define your email 
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=yarmola@princeton.edu

base_dir="/u/yarmola/margulis/margulis-search"
bin_dir="$base_dir/bin"
merge="$base_dir/scripts/merge_trees.py"
data_dir="/scratch/network/yarmola/margulis"
source="output600"
output="validate600"

src1="$data_dir/output600"
out="$data_dir/clean600"

cd $bin_dir

python3 "$merge" -s 1000000 -d 30 -c 1 "$out" "$src1" > "$data_dir/merge.log" 2>&1
