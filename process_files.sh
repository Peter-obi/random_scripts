#!/bin/bash
# Loop through all files ending with colvar.dat in the specified directory
for file in ./*colvar.dat; do
    # Extract the base name without the extension
    base_name=$(basename "$file" .dat)

    # Call the Python script with each file and the corresponding output file name
    python process_file.py "$file" "${base_name}_traj.npy"
done

