#!/bin/bash
#SBATCH -N 1 # node count
#SBATCH --ntasks-per-node=38
#SBATCH -t 239:59:00
#SBATCH --mem=60GB
# sends mail when process begins, and 
# when it ends. Make sure you define your email 
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=yarmola@princeton.edu

base_dir="/u/yarmola/margulis-center/margulis-search"
bin_dir="$base_dir/bin"
words_dir=$base_dir

search="$base_dir/scripts/dotest.py"

data_dir="/scratch/network/yarmola/searchtest"

cd $bin_dir

cat "$data_dir/refine.log" >> "$data_dir/refine.log.all"

python2 "$search" -i 72 -t 12 -r "$bin_dir/search_test" -c 38 "$data_dir/source" "$data_dir/output_live" > "$data_dir/refine.log" 2>&1
