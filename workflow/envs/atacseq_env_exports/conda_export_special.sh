#!/bin/bash


# Path to the input file containing paths
input_file="atacseq_conda_envs.txt"

# Iterate through each line in the file
while IFS= read -r line; do
    # Save the full path into a variable
    full_path="$line"

    # Save the final item in the path to another variable
    # by using the basename command
    final_item=$(basename "$line")

    command_history="conda env export --prefix $full_path --no-build --from-history > ./atacseq_envs/$final_item.fromHistory.yaml"
    command_nobuild="conda env export --prefix $full_path --no-build > ./atacseq_envs/$final_item.noBuild.yaml"
    command_all="conda env export --prefix $full_path > ./atacseq_envs/$final_item.all.yaml"
    
    eval $command_history
    eval $command_nobuild
    eval $command_all
   
done < "$input_file"
     