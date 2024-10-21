import h5py
import argparse

# Predefined order for attributes
attribute_order = ['obs', 'var', 'uns', 'obsm', 'varm', 'layers', 'obsp']

def get_dataset_shape(dataset):
    """
    Safely retrieve the shape of a dataset, including sparse formats.
    """
    # Check if it's a sparse format, typically stored in subgroups like 'data' in CSR or CSC format
    if isinstance(dataset, h5py.Group):
        # Try to handle sparse matrix formats if necessary
        if 'data' in dataset and 'indices' in dataset:
            # Assume this is a sparse matrix (e.g., CSR or CSC)
            shape = dataset.attrs.get('shape')
            if shape is not None:
                return shape
        else:
            raise AttributeError("'X' appears to be a group but doesn't match expected sparse matrix structure.")
    elif isinstance(dataset, h5py.Dataset):
        return dataset.shape
    else:
        raise AttributeError("'X' is neither a recognized dataset nor a sparse group format.")

def display_adata_attributes(file_path):
    # Open the .h5ad file using h5py
    with h5py.File(file_path, "r") as f:
        # Retrieve the shape of the X dataset to get cell and gene count
        if 'X' in f:
            try:
                n_obs, n_vars = get_dataset_shape(f['X'])
                print(f"AnnData object with n_obs × n_vars = {n_obs} × {n_vars}")
            except AttributeError as e:
                print(f"Unable to determine cell and gene count from 'X'. Error: {e}")
        
        # Ensure the attributes are displayed in the specified order
        for attribute in attribute_order:
            if attribute in f:
                # Check if the attribute has columns or subkeys
                if isinstance(f[attribute], h5py.Group):
                    subkeys = list(f[attribute].keys())
                    formatted_columns = ', '.join(f"'{subkey}'" for subkey in subkeys)
                    print(f"{attribute}: {formatted_columns}")
                else:
                    print(f"{attribute}: Not a group, likely an array or other dataset.")
        
        # Display any remaining attributes not listed in the predefined order
        remaining_attributes = set(f.keys()) - set(attribute_order)
        for attribute in remaining_attributes:
            if isinstance(f[attribute], h5py.Group):
                subkeys = list(f[attribute].keys())
                formatted_columns = ', '.join(f"'{subkey}'" for subkey in subkeys)
                print(f"{attribute}: {formatted_columns}")
            else:
                print(f"{attribute}: Not a group, likely an array or other dataset.")

# Set up argument parser
parser = argparse.ArgumentParser(description="Display attributes and columns of an AnnData object (.h5ad file).")
parser.add_argument('file_path', nargs='?', help="The path to the .h5ad file.", default=None)

# Parse arguments
args = parser.parse_args()

# If file path is not provided via command line, ask for input
if args.file_path is None:
    file_path = input("Please enter the path to your .h5ad file: ")
else:
    file_path = args.file_path

# Call the function with the provided or inputted file path
display_adata_attributes(file_path)

