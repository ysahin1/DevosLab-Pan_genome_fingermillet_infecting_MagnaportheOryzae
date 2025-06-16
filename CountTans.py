import os
import pickle
import csv
from collections import Counter

# === Step 1: Setup paths ===
pathtopicklefolder = '../mergedpickle2'           # Folder containing merged .pkl files
asso = '../AssoLocationCounts'                    # Folder to write per-strain interval info
csv_file_path = '../CountsInEachStrain'           # Output folder for summary CSV

# List all .pkl files in the folder
picklelist = []
for file in os.listdir(pathtopicklefolder):
    if file.endswith('.pkl'):
        picklelist.append(os.path.join(pathtopicklefolder, file))

# Dictionary to collect results across all strains
all_strain_result = {}

# === Step 2: Process each .pkl file ===
for picklepath in picklelist:
    # Extract strain name from filename
    picklefile = picklepath.strip().split('/')[-1]
    readname1 = picklefile.strip().split('_')[2]
    strain_name = readname1.strip().split('.')[0]
    file_name = f'{strain_name}'

    # Open output .txt file for writing interval info
    with open(f'{asso}/{file_name}.txt', 'w') as asso1:
        try:
            # Load dictionary from .pkl
            with open(picklepath, 'rb') as pickle_file_1:
                dic = pickle.load(pickle_file_1)

            # Extract coordinate values (column index 1) for difference calculation
            columnC = [dic[id][1] for id in dic]

            # Compute pairwise differences between consecutive positions
            diff = []
            for i, id in enumerate(dic):
                if i == 0:
                    diff.append(columnC[0])
                else:
                    dif = int(columnC[i]) - int(columnC[i - 1])
                    diff.append(dif)

            # Function to find intervals where values in a list remain constant
            def find_intervals(data):
                intervals = []
                start = 0
                for i in range(1, len(data)):
                    if data[i] != data[i - 1]:
                        intervals.append([data[i - 1], start, i - 1])
                        start = i
                intervals.append([data[-1], start, len(data) - 1])
                return intervals

            # Get CIGAR column (index 10) to identify structural states
            columnL = [dic[id][10] for id in dic]
            intervals = find_intervals(columnL)

            # Get chromosome info (column index 0)
            chr_columnB = [dic[id][0] for id in dic]

            # === Step 3: Analyze each interval ===
            for val, start, end in intervals:
                chr = tuple(set(chr_columnB[start:end]))  # Get chromosome(s) involved in this interval
                asso1.write(f'for {val} difference intervals in {chr}; {diff[start:end+1]}\n')

                # Only analyze intervals of sufficient length
                if len(diff[start:end]) >= 3:
                    if strain_name not in all_strain_result:
                        all_strain_result[strain_name] = {}

                    if chr not in all_strain_result[strain_name]:
                        all_strain_result[strain_name][chr] = {'both': 0, 'single': 0}

                    # Count how many differences >20 occur in the interval
                    if len([i for i in diff[start+1:end] if i > 20]) > 0:
                        all_strain_result[strain_name][chr]['both'] += 1
                    else:
                        all_strain_result[strain_name][chr]['single'] += 1

        except Exception as e:
            print(f"Error processing file {file_name}: {e}")

# === Step 4: Export summarized counts for each strain and chromosome ===
with open(f'{csv_file_path}/all_strainAssocount_result.csv', mode='w', newline='') as file:
    # Write header for each chromosome and type
    file.write('strain_name\tChr01_single_end\tChr01_both_end\tChr02_single_end\tChr02_both_end\tChr03_single_end\tChr03_both_end\tChr04_single_end\tChr04_both_end\tChr05_single_end\tChr05_both_end\tChr06_single_end\tChr06_both_end\tChr07_single_end\tChr07_both_end\n')

    for strain_name, chr_dict in all_strain_result.items():
        output_row = [strain_name]
        for i in range(1, 8):  # Loop over chromosomes Chr01 to Chr07
            chr_key = f'Chr{i:02}'
            chr_tuple = (chr_key,)
            if chr_tuple in chr_dict:
                output_row.append(str(chr_dict[chr_tuple]['single']))
                output_row.append(str(chr_dict[chr_tuple]['both']))
            else:
                output_row.append('0')
                output_row.append('0')
        file.write('\t'.join(output_row) + '\n')

