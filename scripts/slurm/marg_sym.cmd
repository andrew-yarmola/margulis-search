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

base_dir="/u/yarmola/margsym/margulis-search"
bin_dir="$base_dir/bin"
words_dir=$base_dir

search="$base_dir/scripts/dosearch.py"
words="$words_dir/words"
impossible="$words_dir/impossible"

data_dir="/scratch/network/yarmola/marg_sym_param"
log_file="refine.log"
output="output_live"

cd $bin_dir

cat "$data_dir/$log_file" >> "$data_dir/${log_file}.all"

python3 "$search" -i 120 -t 12 -r "$bin_dir/refine_marg" -w "$words" -p "$impossible" -c 38 "$data_dir/merged" "$data_dir/$output" > "$data_dir/$log_file" 2>&1
