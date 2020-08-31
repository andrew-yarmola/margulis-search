#!/bin/bash
#SBATCH -N 1 # node count
#SBATCH --ntasks-per-node=10
#SBATCH -t 239:59:00
#SBATCH --mem=60GB
# sends mail when process begins, and 
# when it ends. Make sure you define your email 
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=yarmola@princeton.edu

base_dir="/Users/yarmola/Projects/margulis-search"
bin_dir="$base_dir/bin"
words_dir=$base_dir

search="$base_dir/scripts/dosearch.py"
words="$words_dir/words"
powers="$words_dir/powers"

data_dir="/Users/yarmola/Projects/margulis-search/data"

cd $bin_dir

python2 "$search" -i 42 -t 6 -r "$bin_dir/refine_marg" -w "$words" -p "$powers" -c 6 "$data_dir/source" "$data_dir/output_live" > "$data_dir/refine.log" 2>&1
