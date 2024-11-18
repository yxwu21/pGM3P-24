import argparse
import re
import csv
import numpy as np
from scipy.interpolate import CubicSpline

def extract_data(input_file, output_file, keyword_column_pairs):
    extracted_data = {key: [] for key in keyword_column_pairs.keys()}
    temp0 = None
    
    with open(input_file, 'r') as file:
        lines = file.readlines()

    # Extract temp0
    for line in lines:
        if 'temp0=' in line:
            temp0_match = re.search(r'temp0=(\d+)', line)
            if temp0_match:
                temp0 = float(temp0_match.group(1))
                break

    for key, column in keyword_column_pairs.items():
        pattern = re.compile(rf"{key}\s+")
        for line in reversed(lines):
            if pattern.search(line):
                try:
                    extracted_data[key].append(line.split()[column - 1])
                except IndexError:
                    pass
                if len(extracted_data[key]) == 2:
                    break

    # Print the extracted data
    for key, values in extracted_data.items():
        print(f"{key}: {values}")
    if temp0:
        print(f"Temperature (temp0): {temp0}")

    # Save the extracted data to CSV
    with open(output_file, 'w', newline='') as csvfile:
        fieldnames = ['Property', 'Average', 'STD']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

        writer.writeheader()
        for key, values in extracted_data.items():
            if len(values) == 2:
                writer.writerow({'Property': key, 'Average': values[0], 'STD': values[1]})
            else:
                writer.writerow({'Property': key, 'Average': 'N/A', 'STD': 'N/A'})
        if temp0:
            writer.writerow({'Property': 'temp0', 'Average': temp0, 'STD': 'N/A'})

    print(f"Data extracted and saved to {output_file}")

def calculate_properties(input_file, output_file):
    data = {}
    temp0 = None

    with open(input_file, 'r') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            if row['Property'] == 'temp0':
                temp0 = float(row['Average'])
            else:
                data[row['Property']] = {'avg': float(row['Average']), 'std': float(row['STD'])}

    if not temp0:
        print("Temperature (temp0) not found in the extracted data.")
        return

    density_avg = data['Density']['avg']
    ep_std = data['EPtot']['std']
    vol_avg = data['VOLUME']['avg']
    vol_std = data['VOLUME']['std']

    isobaric_heat_capacity = (ep_std**2) / (1.9872041E-3 * temp0**2) + 1.98
    isothermal_compressibility = (vol_std**2 * 1E4) / (vol_avg * 1.380649 * temp0)

    print(f"Temperature: {temp0}")
    print(f"Density: {density_avg}")
    print(f"Isobaric Heat Capacity: {isobaric_heat_capacity}")
    print(f"Isothermal Compressibility: {isothermal_compressibility}")

    with open(output_file, 'w') as file:
        file.write(f"Temperature: {temp0}\n")
        file.write(f"Density: {density_avg}\n")
        file.write(f"Isobaric Heat Capacity: {isobaric_heat_capacity}\n")
        file.write(f"Isothermal Compressibility: {isothermal_compressibility}\n")

    print(f"Calculated properties saved to {output_file}")

def merge_datasets(input_files, output_file):
    merged_data = []
    for input_file in input_files:
        with open(input_file, 'r') as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                merged_data.append(row)

    # Sort by temperature
    merged_data.sort(key=lambda x: float(x['Average']) if x['Property'] == 'temp0' else float('inf'))

    # Save the merged data to CSV
    with open(output_file, 'w', newline='') as csvfile:
        fieldnames = ['Property', 'Average', 'STD']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for row in merged_data:
            writer.writerow(row)

    print(f"Merged dataset saved to {output_file}")

def calculate_thermal_expansion(input_file, output_file):
    temperatures = []
    densities = []

    with open(input_file, 'r') as csvfile:
        reader = csv.reader(csvfile)
        next(reader)  # Skip header
        for row in reader:
            temperatures.append(float(row[0]))
            densities.append(float(row[1]))

    temperatures = np.array(temperatures)
    densities = np.array(densities)
    
    # Fit density w.r.t. temperature using cubic spline
    cs = CubicSpline(temperatures, densities)
    
    # Calculate the thermal expansion coefficient
    dln_density_dT = cs.derivative()(temperatures) / densities

    with open(output_file, 'w', newline='') as csvfile:
        fieldnames = ['Temperature', 'Thermal Expansion Coefficient']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for temp, exp_coeff in zip(temperatures, dln_density_dT):
            writer.writerow({'Temperature': temp, 'Thermal Expansion Coefficient': exp_coeff})

    print(f"Thermal expansion coefficients saved to {output_file}")

def main():
    parser = argparse.ArgumentParser(description="Process AMBER MD output files.")
    subparsers = parser.add_subparsers(dest="command")

    # Calculate sub-command
    calculate_parser = subparsers.add_parser('calculate', help="Calculate property from input file")
    calculate_parser.add_argument('-i', '--input', required=True, help="Input file to read")
    calculate_parser.add_argument('-o', '--output', required=True, help="Output file to save the result")

    # Extract sub-command
    extract_parser = subparsers.add_parser('extract', help="Extract data from MD output file")
    extract_parser.add_argument('-i', '--input', required=True, help="Input file to read")
    extract_parser.add_argument('-o', '--output', required=True, help="Output file to save the result")

    # Merge sub-command
    merge_parser = subparsers.add_parser('merge', help="Merge datasets from multiple files")
    merge_parser.add_argument('-i', '--input', required=True, nargs='+', help="Input files to read")
    merge_parser.add_argument('-o', '--output', required=True, help="Output file to save the merged result")

    # Special sub-command
    special_parser = subparsers.add_parser('special', help="Calculate special properties")
    special_parser.add_argument('-p', '--property', required=True, help="Property to calculate")
    special_parser.add_argument('-i', '--input', required=True, help="Input file to read")
    special_parser.add_argument('-o', '--output', required=True, help="Output file to save the result")

    args = parser.parse_args()

    keyword_column_pairs = {
        'VOLUME': 9,
        'EPtot': 9,
        'Density': 3
    }

    if args.command == 'calculate':
        calculate_properties(args.input, args.output)
    elif args.command == 'extract':
        extract_data(args.input, args.output, keyword_column_pairs)
    elif args.command == 'merge':
        merge_datasets(args.input, args.output)
    elif args.command == 'special':
        if args.property == 'thermal_expansion':
            calculate_thermal_expansion(args.input, args.output)
        else:
            print("Special property calculation not implemented for this property.")
    else:
        parser.print_help()

if __name__ == "__main__":
    main()

