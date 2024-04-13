import pandas as pd
import numpy as np
import sys

# Get file names from command line arguments
input_file = sys.argv[1]
output_file = sys.argv[2]

# Read the data from the input file
df = pd.read_csv(input_file, sep='\s+', header=None)

# Save the data to a numpy file
np.save(output_file, df)

