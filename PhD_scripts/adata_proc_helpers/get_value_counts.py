import h5py
import pandas as pd
import argparse
import numpy as np

def get_value_counts_from_obs(file_path, obs_column):
    """
    Function to get value counts from a specific column in obs without loading full AnnData.
    """
    with h5py.File(file_path, "r") as f:
        # Load only the obs data
        if 'obs' in f:
            obs_group = f['obs']
            
            # Check if the specified column exists in obs
            if obs_column.encode('utf-8') in obs_group:
                # If the column is a group, handle it as a categorical variable
                if isinstance(obs_group[obs_column.encode('utf-8')], h5py.Group):
                    group = obs_group[obs_column.encode('utf-8')]
                    
                    # Load categories and codes
                    categories = group['categories'][:]
                    codes = group['codes'][:]
                    
                    # Decode byte strings if necessary
                    if isinstance(categories[0], bytes):
                        categories = np.array([x.decode('utf-8') for x in categories])
                    
                    # Map codes to categories
                    obs_data = pd.Series(categories[codes])
                else:
                    # If it's not a group, assume it's a regular dataset
                    obs_data = obs_group[obs_column.encode('utf-8')][:]
                    
                    # Convert byte strings to normal strings if necessary
                    if isinstance(obs_data[0], bytes):
                        obs_data = pd.Series([x.decode('utf-8') for x in obs_data])
                    else:
                        obs_data = pd.Series(obs_data)
                
                # Get value counts
                value_counts = obs_data.value_counts()
                
                # Print value counts
                print(f"Value counts for '{obs_column}':")
                print(value_counts)
            else:
                print(f"Column '{obs_column}' not found in obs.")
        else:
            print("'obs' group not found in the file.")

# Set up argument parser
parser = argparse.ArgumentParser(description="Get value counts from a specific column in obs from an AnnData object (.h5ad file) without loading full data.")
parser.add_argument('file_path', help="The path to the input .h5ad file.")
parser.add_argument('obs_column', help="The column in obs to get value counts from (e.g., 'study').")

# Parse arguments
args = parser.parse_args()

# Call the function with the provided arguments
get_value_counts_from_obs(args.file_path, args.obs_column)

