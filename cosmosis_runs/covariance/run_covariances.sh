#!/bin/bash

script_file="run_covariance.py"

for i in {0..3}; do
    echo "Bin$i"
    while read leaf_id; do
        output_dir="/global/cfs/cdirs/des/elisa/ShearSplits_data/decisiontree/cosmosis_runs/bin${i}/leaf_${leaf_id}"
        if (( leaf_id % 100 == 0 )); then
            echo "  Running leaf $leaf_id"
        fi
        python3 "$script_file" "$output_dir" &
        
        # Optional: Limit number of parallel jobs
        while [ "$(jobs -r | wc -l)" -ge 4 ]; do
            sleep 1
        done
    done < "leaf_ids_bin${i}.txt"
done

# Wait for all background jobs to finish
wait
echo "All jobs completed."