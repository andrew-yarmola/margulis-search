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

#base_dir="/Users/yarmola/Projects/margulis-search"
base_dir="/u/yarmola/margulis-search"
bin_dir="$base_dir/bin"
words_dir=$base_dir

search="$base_dir/scripts/dosearch.py"
words="$words_dir/words"
powers="$words_dir/powers"

#data_dir="/Users/yarmola/Projects/margulis-search/data"
data_dir="/scratch/network/yarmola/margulis"

cd $bin_dir

cat "$data_dir/refine.log" >> "$data_dir/refine.log.all"

python2 "$search" -i 22 -t 6 -r "$bin_dir/refine_marg" -w "$words" -p "$powers" -c 38 "$data_dir/source" "$data_dir/output_live" > "$data_dir/refine.log" 2>&1
