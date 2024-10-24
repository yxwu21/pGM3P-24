import pandas as pd
import glob
import os


def merge_data(existing_data, new_data):
    """
    Merges new data into the existing structured data dictionary.

    Parameters:
    existing_data (dict): The existing structured data.
    new_data (dict): The new data to be merged.

    Returns:
    dict: The updated existing data with the new data merged in.
    """
    for prefix in new_data:
        if prefix not in existing_data:
            existing_data[prefix] = new_data[prefix]
        else:
            for sub_key in new_data[prefix]:
                # Merge sub_key data
                if sub_key not in existing_data[prefix]:
                    existing_data[prefix][sub_key] = new_data[prefix][sub_key]
                else:
                    # Merge step data
                    existing_data[prefix][sub_key].update(new_data[prefix][sub_key])
    return existing_data


def structure_data(df):
    """
    Function to structure the DataFrame into a nested dictionary format.
    The first level of keys will be based on unique prefixes in the column names.
    The second level of keys will be the sub-keys extracted from the column names.
    The values will be dictionaries mapping 'Step' to the corresponding data.

    Parameters:
    df (pandas.DataFrame): The DataFrame to be structured.

    Returns:
    dict: A nested dictionary with the structured data.
    """
    structured_data = {}

    # Iterate over each column (excluding 'Step')
    for col in df.columns[1:]:
        # Splitting the column name to extract the prefix and the sub-key
        parts = col.split(" - ")
        if len(parts) > 1:
            prefix, sub_key = parts
            sub_key = sub_key.split("/")[
                -1
            ]  # Further splitting to get the exact sub-key

            # Initialize the nested dictionary if not already done
            if prefix not in structured_data:
                structured_data[prefix] = {}

            # Populate the data
            structured_data[prefix][sub_key] = df.set_index("Step")[col].to_dict()

    return structured_data


def process_all_csv_files(directory_path):
    """
    Processes all CSV files in the given directory.

    Parameters:
    directory_path (str): Path to the directory containing CSV files.

    Returns:
    dict: A dictionary with structured data from all CSV files.
    """
    all_data = {}
    for file_path in glob.glob(os.path.join(directory_path, "*.csv")):
        df = pd.read_csv(file_path)
        structured_data = structure_data(df)
        all_data = merge_data(all_data, structured_data)
    return all_data


# Example usage
if __name__ == "__main__":
    directory_path = "/home/yxwu/pGM_water_model/wandb_data/vital-cosmos-120"
    wandb_data = process_all_csv_files(directory_path)

    # Example: Accessing data
    step = 10
    scale_p_value = (
        wandb_data.get("vital-cosmos-120", {})
        .get("config.scale_p", {})
        .get(step, "No data for this step")
    )
    # print(f"Scale_p value at Step {step}: {scale_p_value}")
    print(wandb_data.get("vital-cosmos-120", {}).keys())
