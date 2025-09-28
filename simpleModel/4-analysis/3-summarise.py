import os
import csv
import re
import statistics

def parse_filename(filename):
    """Extract r, minD, and minVar from filename"""
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
    """Read a stats CSV file and return metrics as dict of lists"""
    metrics = {
        'TPR': [], 'FPR': [], 'FNR': [], 'TNR': [],
        'Accuracy': [], 'Precision': [], 'F1Score': []
    }
    
    try:
        with open(filepath, 'r') as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                try:
                    for key in metrics.keys():
                        metrics[key].append(float(row[key]))
                except (ValueError, KeyError) as e:
                    print(f"Warning: Error parsing row in {filepath}: {e}")
                    continue
                    
    except Exception as e:
        print(f"Error reading {filepath}: {e}")
        return {k: [] for k in metrics}
    
    return metrics

def calculate_means(values):
    """Calculate mean, handling empty lists"""
    if not values:
        return 0.0
    return statistics.mean(values)

def main():
    stats_folder = "./stats"
    output_folder = "./summary"
    output_file = os.path.join(output_folder, "summary.csv")

    os.makedirs(output_folder, exist_ok=True)
    
    if not os.path.exists(stats_folder):
        print(f"Stats folder {stats_folder} does not exist!")
        return
    
    summary_data = []
    
    # Process all CSV files in stats folder
    for filename in os.listdir(stats_folder):
        if not filename.startswith("stats_r") or not filename.endswith(".csv"):
            continue
            
        r, minD, minVar = parse_filename(filename)
        if r is None:
            print(f"Could not parse filename: {filename}")
            continue
        
        filepath = os.path.join(stats_folder, filename)
        metrics = read_stats_file(filepath)
        
        if not metrics['TPR']:
            print(f"No valid data found in {filename}")
            continue
        
        # Calculate mean for each metric
        mean_metrics = {f"Mean_{k}": calculate_means(v) for k, v in metrics.items()}
        
        # Store summary data
        row = {
            'r_value': r,
            'minD_value': minD,
            'minVar_value': minVar
        }
        row.update(mean_metrics)
        summary_data.append(row)
        
        print(f"Processed {filename}: r={r}, minD={minD}, minVar={minVar}")
    
    # Sort by r value
    summary_data.sort(key=lambda x: x['r_value'])
    
    # Write summary CSV
    if summary_data:
        fieldnames = ['r_value', 'minD_value', 'minVar_value'] + list(summary_data[0].keys())[3:]
        
        with open(output_file, 'w', newline='') as csvfile:
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
