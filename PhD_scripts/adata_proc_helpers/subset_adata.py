import anndata as ad
import argparse

def subset_adata(file_path, obs_column, obs_value, output_path):
    """
    Function to subset an AnnData object based on a specific value in a given obs column.
    """
    # Load the AnnData object
    adata = ad.read_h5ad(file_path)
    
    # Check if the obs_column exists in the AnnData object
    if obs_column not in adata.obs.columns:
        print(f"Error: '{obs_column}' not found in obs columns.")
        return
    
    # Subset the AnnData object to include only rows where obs_column matches obs_value
    subset_adata = adata[adata.obs[obs_column] == obs_value, :]
    
    # Save the subset AnnData object to a new file
    subset_adata.write(output_path)
    print(f"Subset AnnData object saved to {output_path}")

# Set up argument parser
parser = argparse.ArgumentParser(description="Subset an AnnData object based on a specific value in a given obs column.")
parser.add_argument('file_path', help="The path to the input .h5ad file.")
parser.add_argument('obs_column', help="The column in obs to filter by (e.g., 'study').")
parser.add_argument('obs_value', help="The specific value in the obs column to filter by (e.g., 'Kaminski_2020').")
parser.add_argument('output_path', help="The path to save the output subset .h5ad file.")

# Parse arguments
args = parser.parse_args()

# Call the function with the provided arguments
subset_adata(args.file_path, args.obs_column, args.obs_value, args.output_path)

