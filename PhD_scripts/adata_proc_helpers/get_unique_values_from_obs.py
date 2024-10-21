import h5py
import argparse
import numpy as np

def get_unique_values_from_obs(file_path, obs_column):
    """
    Function to get unique values from a specific column in obs without fully loading the file.
    """
    with h5py.File(file_path, "r") as f:
        # Check if 'obs' exists in the file
        if 'obs' in f:
            obs_group = f['obs']
            #print(f"Available keys in 'obs': {list(obs_group.keys())}")
            
            # Check if the specified column exists in obs
            if obs_column in obs_group:
                # Check if it's a categorical variable
                if isinstance(obs_group[obs_column], h5py.Group):
                    group = obs_group[obs_column]
                    
                    # Access 'categories' and 'codes'
                    categories = group['categories'][:]
                    codes = group['codes'][:]
                    
                    # Convert byte strings to normal strings if necessary
                    if isinstance(categories[0], bytes):
                        categories = np.array([x.decode('utf-8') for x in categories])
                    
                    # Get unique categories based on codes
                    unique_values = np.unique(categories[codes])
                    print(f"obs['{obs_column}']: {unique_values}")
                else:
                    print(f"Column '{obs_column}' is not stored as a categorical group.")
            else:
                print(f"Column '{obs_column}' not found in obs.")
        else:
            print("'obs' group not found in the file.")

# Set up argument parser
parser = argparse.ArgumentParser(description="Get unique values from a specific column in obs from an AnnData object (.h5ad file).")
parser.add_argument('file_path', help="The path to the .h5ad file.")
parser.add_argument('obs_column', help="The column in obs to extract unique values from (e.g., 'study').")

# Parse arguments
args = parser.parse_args()

# Call the function with the provided arguments
get_unique_values_from_obs(args.file_path, args.obs_column)

