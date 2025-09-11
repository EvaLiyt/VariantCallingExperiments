import os
import csv
import re
import statistics

def parse_filename(filename):
    """Extract r, minD, and minVar from filename"""
    # Pattern: stats_r<number>_minD<number>_minVar<number>.csv
    pattern = r"stats_r(\d+)_minD(\d+)_minVar(\d+)\.csv"
    match = re.match(pattern, filename)
    
    if match:
        r = int(match.group(1))
        minD = int(match.group(2))
        minVar = int(match.group(3))
        return r, minD, minVar
    else:
        return None, None, None

def read_stats_file(filepath):
    """Read a stats CSV file and return TPR, FPR, Precision values"""
    tpr_values = []
    fpr_values = []
    precision_values = []
    
    try:
        with open(filepath, 'r') as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                # Convert string values to float, handle potential errors
                try:
                    tpr = float(row['TPR'])
                    fpr = float(row['FPR'])
                    precision = float(row['Precision'])
                    
                    tpr_values.append(tpr)
                    fpr_values.append(fpr)
                    precision_values.append(precision)
                except (ValueError, KeyError) as e:
                    print(f"Warning: Error parsing row in {filepath}: {e}")
                    continue
                    
    except Exception as e:
        print(f"Error reading {filepath}: {e}")
        return [], [], []
    
    return tpr_values, fpr_values, precision_values

def calculate_means(values):
    """Calculate mean, handling empty lists"""
    if not values:
        return 0.0
    return statistics.mean(values)

def main():
    stats_folder = "./stats"
    output_folder = "./summary"
    output_file = os.path.join(output_folder, "summary.csv")
    output_folder = "./summary"
    output_file = os.path.join(output_folder, "summary.csv")

    # Create folder if it doesn’t exist
    os.makedirs(output_folder, exist_ok=True)
    
    if not os.path.exists(stats_folder):
        print(f"Stats folder {stats_folder} does not exist!")
        return
    
    summary_data = []
    
    # Process all CSV files in the stats folder
    for filename in os.listdir(stats_folder):
        if not filename.startswith("stats_r") or not filename.endswith(".csv"):
            continue
            
        # Parse filename to extract parameters
        r, minD, minVar = parse_filename(filename)
        if r is None:
            print(f"Could not parse filename: {filename}")
            continue
        
        # Read the stats file
        filepath = os.path.join(stats_folder, filename)
        tpr_values, fpr_values, precision_values = read_stats_file(filepath)
        
        if not tpr_values:
            print(f"No valid data found in {filename}")
            continue
        
        # Calculate means
        mean_tpr = calculate_means(tpr_values)
        mean_fpr = calculate_means(fpr_values)
        mean_precision = calculate_means(precision_values)
        
        # Store summary data
        summary_data.append({
            'r_value': r,
            'minD_value': minD,
            'minVar_value': minVar,
            'Mean_TPR': mean_tpr,
            'Mean_FPR': mean_fpr,
            'Mean_Precision': mean_precision
        })
        
        print(f"Processed {filename}: r={r}, minD={minD}, minVar={minVar}")
    
    # Sort by r value for consistent output
    summary_data.sort(key=lambda x: x['r_value'])
    
    # Write summary CSV
    if summary_data:
        with open(output_file, 'w', newline='') as csvfile:
            fieldnames = ['r_value', 'minD_value', 'minVar_value', 'Mean_TPR', 'Mean_FPR', 'Mean_Precision']
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            
            writer.writeheader()
            for row in summary_data:
                writer.writerow(row)
        
        print(f"\nSummary written to {output_file}")
        print(f"Processed {len(summary_data)} files successfully")
    
    else:
        print("No valid stats files found to process!")

if __name__ == "__main__":
    main()
