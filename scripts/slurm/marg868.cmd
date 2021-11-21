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
bad_relators="$words_dir/bad_relators"

data_dir="/scratch/network/yarmola/marg_sym_param"
log_file="refine868.log"
source="source"
output="output868"

cosh_marg="1.4"
sinh_rad="2.2447"

cd $bin_dir

cat "$data_dir/$log_file" >> "$data_dir/${log_file}.all"

python3 "$search" \
  -i 18 \
  -t 6 \
  -r "$bin_dir/refine_marg" \
  -w "$words" \
  -p "$impossible" \
  -b "$bad_relators" \
  -M "$cosh_marg" \
  -R "$sinh_rad" \
  -c 38 \
  "$data_dir/$source" "$data_dir/$output" > "$data_dir/$log_file" 2>&1
