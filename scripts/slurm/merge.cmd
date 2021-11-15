#!/bin/bash
#SBATCH -N 1 # node count
#SBATCH --ntasks-per-node=38
#SBATCH -t 239:59:00
#SBATCH --mem=120GB
# sends mail when process begins, and 
# when it ends. Make sure you define your email 
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=yarmola@princeton.edu

base_dir="/u/yarmola/margsym/margulis-search/"
bin_dir="$base_dir/bin"
merge="$base_dir/scripts/merge_trees.py"
data_dir="/scratch/network/yarmola/marg_sym_param"

src1="/scratch/network/yarmola/marg_sym_param/output_live"

cd $bin_dir

python3 "$merge" -s 1000000 -d 30 -c 38 "$data_dir/merge_live" "$src1" > "$data_dir/merge.log" 2>&1
