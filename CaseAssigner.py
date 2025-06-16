"""
Genomic Breakpoint Analysis Script
==================================

This script processes genomic sequencing data to identify and analyze breakpoints
in paired-end reads, particularly focusing on transposable element insertions.

The pipeline consists of two main stages:
1. Initial processing: Merge R1/R2 reads and perform initial breakpoint analysis
2. Enhanced processing: Add distance information and refine breakpoint calls

Author: [Your Name]
Date: [Date]
"""

import pickle
import csv
import os
import re

# =============================================================================
# CONFIGURATION AND SETUP
# =============================================================================

# Define input and output paths
pathtopicklefolder = '../filteredpickle'  # Directory containing filtered pickle files
pathtoOutputPickle = '../mergedpickle'    # Output directory for merged pickles
csv_file_path = '../infoResult'           # Output directory for CSV results
distance = '../distance/distance.txt'    # Distance reference file

# Build list of pickle files to process
picklelist = []
for file in os.listdir(pathtopicklefolder):
    if file.endswith('.pkl'):
        picklelist.append(os.path.join(pathtopicklefolder, file))

# Load distance data for genomic coordinates
with open(distance) as fff:
    distance_data = [line.strip('\n').strip().split('\t') for line in fff]

# =============================================================================
# HELPER FUNCTIONS FOR GENOMIC ANALYSIS
# =============================================================================

def sum_m_and_d(categorized_dict):
    """
    Calculate sum of M and D values for genomic alignment scoring.
    
    For M values: Sum all except the maximum (removes largest alignment)
    For D values: Sum all values (deletions)
    
    Args:
        categorized_dict: Dictionary with 'M', 'I', 'D' keys containing lists of integers
        
    Returns:
        int: Total score combining M and D values
    """
    if categorized_dict['M']:
        max_m = max(categorized_dict['M'])  # Find the largest alignment
        sum_m = sum(categorized_dict['M']) - max_m  # Subtract largest from total
    else:
        sum_m = 0

    # Handle 'D': Sum all deletion lengths
    sum_d = sum(categorized_dict['D'])
    total = int(sum_m + sum_d)
    return total

def categorize_numbers(s):
    """
    Parse CIGAR-like string to categorize alignment operations.
    
    Extracts numbers preceding M (match), I (insertion), D (deletion) characters
    and organizes them into separate lists.
    
    Args:
        s: String containing pattern like "100M50I25D"
        
    Returns:
        dict: Dictionary with 'M', 'I', 'D' keys containing lists of integers
    """
    categorized_numbers = {'M': [], 'I': [], 'D': []}
    
    # Find all number-letter pairs in the string
    for match in re.finditer(r'(\d+)([MID])', s):
        number = int(match.group(1))  # Extract the number
        letter = match.group(2)       # Extract the operation (M/I/D)
        categorized_numbers[letter].append(number)
    
    return categorized_numbers

def extract_numbers(part):
    """
    Extract all numeric values from a string.
    
    Args:
        part: String that may contain numbers
        
    Returns:
        list: List of integers found in the string
    """
    return [int(num) for num in re.findall(r'\d+', part)]

def compare_numbers(s):
    """
    Analyze alignment string to determine breakpoint characteristics.
    
    Splits alignment string by M/I/D operations and determines whether
    the largest alignment is on the left or right side, then calculates
    appropriate scoring based on position.
    
    Args:
        s: Alignment string (CIGAR-like format)
        
    Returns:
        list: List containing tuple with (original_string, side, max_value, score)
    """
    results = []
    
    # Split the string by alignment operations
    parts = re.split(r'[MID]', s)
    
    # Separate numbers into left and right groups
    left_numbers = []
    right_numbers = []
    for i, part in enumerate(parts):
        if i == 0:
            left_numbers.extend(extract_numbers(part))
        else:
            numbers = extract_numbers(part)
            if len(left_numbers) > 0:
                right_numbers.extend(numbers)
            else:
                left_numbers.extend(numbers)
    
    # Get categorized alignment operations
    categorized_numbers = categorize_numbers(s)

    # Find maximum values on each side
    max_left = max(left_numbers) if left_numbers else float('-inf')
    max_right = max(right_numbers) if right_numbers else float('-inf')
    
    # Determine scoring based on which side has larger alignment
    if max_left > max_right:  # Left side dominant
        picknumber = max(categorized_numbers['M'])  # Use largest M value
        result_tuple = (s, 'left', max_left, picknumber)
        results.append(result_tuple)
    else:  # Right side dominant
        picknumber = sum_m_and_d(categorized_numbers)  # Use combined M+D score
        result_tuple = (s, 'right', max_right, picknumber)
        results.append(result_tuple)
    
    return results

# =============================================================================
# STAGE 1: INITIAL BREAKPOINT PROCESSING
# =============================================================================

print("Starting Stage 1: Initial breakpoint processing...")

for picklepath in picklelist:
    # Extract file information from path
    picklefile = picklepath.strip().split('/')[-1]
    readname1 = picklefile.strip().split('_')[2]
    readname = readname1.strip().split('.')[0]
    strain_name = picklefile.strip().split('_')[1]
    file_name = f'{strain_name}_{readname}'

    # Define output paths
    file_path = f'{csv_file_path}/bp_result_{strain_name}_assgined.csv'
    final_pickle_file = f'{pathtoOutputPickle}/mr_{strain_name}.pkl'

    try:
        # Load R1 reads
        with open(picklepath, 'rb') as pickle_file_1:
            read1 = pickle.load(pickle_file_1)
        
        # Only process R1 files (paired-end sequencing)
        if readname != 'R1':
            continue
        else:
            # Mark all R1 reads with identifier '1'
            for id in read1:
                read1[id].insert(0, '1')
            
            # Try to load corresponding R2 file
            check_file = f'{pathtopicklefolder}/info2_{strain_name}_R2.pkl'
            if os.path.exists(check_file):
                with open(check_file, 'rb') as pickle_file_2:
                    read2 = pickle.load(pickle_file_2)
                    # Mark all R2 reads with identifier '2'
                    for id in read2:
                        read2[id].insert(0, '2')
            else:
                continue

        # Merge R1 and R2 reads into single dataset
        mr = {**read1, **read2}

        # =============================================================================
        # BREAKPOINT ANALYSIS FOR EACH READ PAIR
        # =============================================================================
        
        for id in mr:
            # Extract column data (genomic coordinates and alignment info)
            columnD = int(mr[id][2])    # Start position 1
            columnE = int(mr[id][3])    # Start position 2  
            columnF = str(mr[id][4])    # Quality flag
            columnG = int(mr[id][5])    # Reference start
            columnH = int(mr[id][6])    # Reference end 1
            columnI = int(mr[id][7])    # Reference end 2
            columnJ = int(mr[id][8])    # Length
            columnL = int(mr[id][10])   # Target start
            columnM = str(mr[id][11])   # CIGAR string 1
            columnN = str(mr[id][12])   # Direction 1
            columnO = str(mr[id][13])   # CIGAR string 2
            columnP = str(mr[id][14])   # Direction 2

            # CASE 1: Potential deletion/inversion (E < D and quality OK)
            if columnE < columnD and columnF == 'Yes':
                if columnN == 'right':
                    # Check for match-only alignment (no insertions/deletions)
                    if 'M' in columnO and 'I' not in columnO or 'D' not in columnO:
                        x1 = compare_numbers(columnM)[0][3]
                        y1 = columnG + x1
                        # Check if breakpoint falls within expected range
                        if columnL - 15 < y1 < columnL + columnJ + 15:
                            mr[id].insert(0, '1/4rev')  # Classify as type 1/4 reverse
                            mr[id].insert(0, y1)
                        else:
                            mr[id].insert(0, 'None')
                            mr[id].insert(0, y1)
                    else:
                        mr[id].insert(0, 'None')
                        mr[id].insert(0, 'None')
                        
                elif columnP == 'left':  # Case 4: Different orientation
                    if 'M' in columnM and 'I' not in columnM or 'D' not in columnM:
                        if len(mr[id]) > 16:
                            columnR = int(mr[id][16])  # Additional coordinate
                            x3 = compare_numbers(columnO)[0][3]
                            y3 = columnR + x3
                            if columnI - columnH - 15 < y3 < columnI + 15:
                                mr[id].insert(0, '1rev/4')
                                mr[id].insert(0, y3)
                            else:
                                mr[id].insert(0, 'None')
                                mr[id].insert(0, y3)
                        else:
                            mr[id].insert(0, 'None')
                            mr[id].insert(0, 'None')
                    else:
                        mr[id].insert(0, 'None')
                        mr[id].insert(0, 'None')
                else:
                    mr[id].insert(0, 'None')
                    mr[id].insert(0, 'None')
            
            # CASE 2: Mate pair maps to different location (E > D and quality OK)
            elif columnE > columnD and columnF == 'Yes':
                if columnN == 'right':
                    if columnO == 'None':  # No second alignment
                        if columnL - 15 < columnG < columnL + columnJ + 15:
                            mr[id].insert(0, '2rev/3')
                            mr[id].insert(0, columnG)
                        else:
                            mr[id].insert(0, 'None')
                            mr[id].insert(0, 'None')
                    else:
                        mr[id].insert(0, 'None')
                        mr[id].insert(0, 'None')
                        
                elif columnP == 'right':  # Case 3: Mate pair doesn't map
                    if columnM == '*':  # No first alignment
                        if len(mr[id]) > 16:
                            columnR = int(mr[id][16])
                            x2 = compare_numbers(columnO)[0][3]
                            y2 = columnR + x2
                            if columnL - 15 < y2 < columnL + columnJ + 15:
                                mr[id].insert(0, '2rev/3')
                                mr[id].insert(0, y2)
                            else:
                                mr[id].insert(0, 'None')
                                mr[id].insert(0, y2)
                        else:
                            mr[id].insert(0, 'None')
                            mr[id].insert(0, 'None')
                    else:
                        mr[id].insert(0, 'None')
                        mr[id].insert(0, 'None')
                        
                elif columnN == 'left':  # Case 5: Left orientation
                    if columnO == 'None':
                        y4 = columnG + columnD
                        if columnI - columnH - 15 < y4 < columnI + 15:
                            mr[id].insert(0, '2/3rev')
                            mr[id].insert(0, y4)
                        else:
                            mr[id].insert(0, 'None')
                            mr[id].insert(0, y4)
                    else:
                        mr[id].insert(0, 'None')
                        mr[id].insert(0, 'None')
                else:
                    mr[id].insert(0, 'None')
                    mr[id].insert(0, 'None')
            else:
                # No valid breakpoint detected
                mr[id].insert(0, 'None')
                mr[id].insert(0, 'None')

        # Save processed data to pickle file
        try:
            with open(final_pickle_file, 'wb') as outpicklefile:
                pickle.dump(mr, outpicklefile)
        except Exception as e:
            print(f"Error creating pickle file mr_{strain_name}.pkl: {e}")

        # Export to CSV for human readability
        try:
            with open(file_path, 'w', newline='') as csvfile:
                # Create dynamic column headers based on maximum data length
                fieldnames = ['Key'] + [f'Value_{i+1}' for i in range(max(len(values) for values in mr.values()))]
                csv_writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
                csv_writer.writeheader()
                
                # Write each read's data as a row
                for key, values in mr.items():
                    row_data = {'Key': key}
                    for i, value in enumerate(values, start=1):
                        row_data[f'Value_{i}'] = value
                    csv_writer.writerow(row_data)
        except Exception as e:
            print(f"Error creating CSV file: {file_path}: {e}")

    except Exception as e:
        print(f"Error processing file {picklepath}: {e}")

# =============================================================================
# STAGE 2: ENHANCED PROCESSING WITH DISTANCE INFORMATION
# =============================================================================

print("Starting Stage 2: Enhanced processing with distance information...")

# Setup for second stage processing
pathtoOutputPickle_2 = '../mergedpickle2'
csv_file_path_2 = '../infoResult2'
pathtopicklefolder = '../mergedpickle'  # Input from stage 1 output

# Get list of merged pickle files from stage 1
picklelist = []
for file in os.listdir(pathtopicklefolder):
    if file.endswith('.pkl'):
        picklelist.append(os.path.join(pathtopicklefolder, file))

for picklepath in picklelist:
    # Extract strain information
    picklefile = picklepath.strip().split('/')[-1]
    readname1 = picklefile.strip().split('_')[1]
    strain_name = readname1.strip().split('.')[0]
    file_name = f'mr_version2_{strain_name}'

    # Define output paths for stage 2
    file_path = f'{csv_file_path_2}/bp_result_{strain_name}_Assigned.csv'
    final_pickle_file = f'{pathtoOutputPickle_2}/{file_name}.pkl'
    
    try:
        # Load processed data from stage 1
        with open(picklepath, 'rb') as pickle_file_2:
            mr_assigned = pickle.load(pickle_file_2)

        # =============================================================================
        # DISTANCE INFORMATION INTEGRATION
        # =============================================================================
        
        # Ensure all records have sufficient columns (pad if necessary)
        for id in mr_assigned:
            lenlist = len(mr_assigned[id])
            if lenlist <= 18:
                # Add placeholder values for missing columns
                mr_assigned[id].insert(18, 'None')
                mr_assigned[id].insert(19, 'None')
                mr_assigned[id].insert(20, 'None')
                mr_assigned[id].insert(21, 'None')
                mr_assigned[id].insert(22, 'None')
                mr_assigned[id].insert(23, 'None')
            else:
                # Extract breakpoint coordinates for distance lookup
                bpstart = str(mr_assigned[id][19])  # Breakpoint start
                bpchr = str(mr_assigned[id][20])    # Chromosome
                
                # Find matching distance information
                for line in distance_data:
                    dischr = str(line[0])  # Distance file chromosome
                    bp1 = line[1]          # Distance file position
                    dis = line[2]          # Distance value
                    
                    if bpstart == bp1 and bpchr == dischr:
                        mr_assigned[id].insert(19, dis)  # Insert distance
                        break

        # Add distance information for end coordinates
        for id in mr_assigned:
            lenlist = len(mr_assigned[id])
            if lenlist > 18:
                bpend = str(mr_assigned[id][22])   # Breakpoint end
                bpchr = str(mr_assigned[id][21])   # Chromosome
                
                # Find distance to next breakpoint
                for i in range(len(distance_data)):
                    dischr = str(distance_data[i][0])
                    bp1 = distance_data[i][1]
                    dis = distance_data[i][2]
                    
                    if bpend == bp1 and bpchr == dischr:
                        # Insert distance from next entry
                        mr_assigned[id].insert(21, distance_data[i+1][2])
                        break

        # =============================================================================
        # REFINED BREAKPOINT CLASSIFICATION
        # =============================================================================
        
        for id in mr_assigned:
            # Extract all relevant columns for analysis
            columnD = int(mr_assigned[id][4])   # Position 1
            columnE = int(mr_assigned[id][5])   # Position 2
            columnG = int(mr_assigned[id][7])   # Reference start
            columnH = int(mr_assigned[id][8])   # Reference end 1
            columnI = int(mr_assigned[id][9])   # Reference end 2
            columnJ = int(mr_assigned[id][10])  # Length
            columnK = str(mr_assigned[id][11])  # Additional info
            columnL = int(mr_assigned[id][12])  # Target start
            columnM = str(mr_assigned[id][13])  # CIGAR 1
            columnN = str(mr_assigned[id][14])  # Direction 1
            columnO = str(mr_assigned[id][15])  # CIGAR 2
            columnP = str(mr_assigned[id][16])  # Direction 2
            columnR = mr_assigned[id][18]       # Additional coordinate
            columnS = mr_assigned[id][19]       # Distance 1
            columnT = mr_assigned[id][20]       # Coordinate
            columnU = mr_assigned[id][21]       # Distance 2
            columnV = mr_assigned[id][22]       # Chromosome info
            columnW = mr_assigned[id][23]       # End coordinate

            # REFINED CASE 1: Deletion/inversion analysis
            if columnE < columnD:
                if columnN == 'right' and columnP == 'None':
                    x1 = compare_numbers(columnM)[0][3]
                    y1 = columnG + x1
                    if columnL - 15 < y1 < columnL + columnJ + 15:
                        mr_assigned[id][1] = '1/4rev'
                        mr_assigned[id][0] = y1
                        mr_assigned[id].insert(0, columnK)
                    else:
                        mr_assigned[id][1] = 'None'
                        mr_assigned[id][0] = y1
                        mr_assigned[id].insert(0, columnK)

                # REFINED CASE 2: Left orientation with distance check
                elif columnP == 'left' and columnN == 'None':
                    x2 = compare_numbers(columnO)[0][3]
                    y2 = int(columnR) + x2
                    if int(columnT) - int(columnS) - 15 < y2 < int(columnT) + 15:
                        mr_assigned[id][1] = '1/4rev'
                        mr_assigned[id][0] = y2
                        mr_assigned[id].insert(0, columnV)
                    else:
                        mr_assigned[id][1] = 'None'
                        mr_assigned[id][0] = y2
                        mr_assigned[id].insert(0, columnV)

                # NEW CASE: Right orientation analysis
                elif columnP == 'right' and columnN == 'None':
                    x21 = compare_numbers(columnO)[0][3]
                    y21 = int(columnR) + x21
                    if int(columnW) - 15 < y21 < int(columnW) + int(columnU) + 15:
                        mr_assigned[id][1] = '1/4rev'
                        mr_assigned[id][0] = y21
                        mr_assigned[id].insert(0, columnV)
                    else:
                        mr_assigned[id][1] = 'None'
                        mr_assigned[id][0] = y21
                        mr_assigned[id].insert(0, columnV)

                # NEW CASE: Left orientation variant
                elif columnN == 'left' and columnP == 'None':
                    x51 = compare_numbers(columnM)[0][3]
                    y51 = columnG + x51
                    if columnI - columnH - 15 < y51 < columnI + 15:
                        mr_assigned[id][1] = '1rev/4'
                        mr_assigned[id][0] = y51
                        mr_assigned[id].insert(0, columnK)
                    else:
                        mr_assigned[id][1] = 'None'
                        mr_assigned[id][0] = y51
                        mr_assigned[id].insert(0, columnK)
                else:
                    mr_assigned[id][1] = 'None'
                    mr_assigned[id][0] = 'None'
                    mr_assigned[id].insert(0, 'None')

            # REFINED CASE 3: Transposable element extension analysis
            elif columnE > columnD:
                # Case 3: Right orientation
                if columnN == 'right' and columnP == 'None':
                    x3 = compare_numbers(columnM)[0][3]
                    y3 = columnG + x3
                    if columnL - 15 < y3 < columnL + columnJ + 15:
                        mr_assigned[id][1] = '2rev/3'
                        mr_assigned[id][0] = y3
                        mr_assigned[id].insert(0, columnK)
                    else:
                        mr_assigned[id][1] = 'None'
                        mr_assigned[id][0] = y3
                        mr_assigned[id].insert(0, columnK)

                # Case 4: Right orientation with distance
                elif columnP == 'right' and columnN == 'None':
                    x4 = compare_numbers(columnO)[0][3]
                    y4 = int(columnR) + x4
                    if int(columnW) - 15 < y4 < int(columnW) + int(columnU) + 15:
                        mr_assigned[id][1] = '2rev/3'
                        mr_assigned[id][0] = y4
                        mr_assigned[id].insert(0, columnV)
                    else:
                        mr_assigned[id][1] = 'None'
                        mr_assigned[id][0] = y4
                        mr_assigned[id].insert(0, columnV)

                # Case 5: Left orientation
                elif columnN == 'left' and columnP == 'None':
                    x5 = compare_numbers(columnM)[0][3]
                    y5 = columnG + x5
                    if columnI - columnH - 15 < y5 < columnI + 15:
                        mr_assigned[id][1] = '2/3rev'
                        mr_assigned[id][0] = y5
                        mr_assigned[id].insert(0, columnK)
                    else:
                        mr_assigned[id][1] = 'None'
                        mr_assigned[id][0] = y5
                        mr_assigned[id].insert(0, columnK)

                # Case 6: Left orientation with distance
                elif columnP == 'left' and columnN == 'None':
                    x6 = compare_numbers(columnO)[0][3]
                    y6 = int(columnR) + x6
                    if int(columnT) - int(columnS) - 15 < y6 < int(columnT) + 15:
                        mr_assigned[id][1] = '2/3rev'
                        mr_assigned[id][0] = y6
                        mr_assigned[id].insert(0, columnV)
                    else:
                        mr_assigned[id][1] = 'None'
                        mr_assigned[id][0] = y6
                        mr_assigned[id].insert(0, columnV)
                else:
                    mr_assigned[id][1] = 'None'
                    mr_assigned[id][0] = 'None'
                    mr_assigned[id].insert(0, 'None')
            else:
                mr_assigned[id][1] = 'None'
                mr_assigned[id][0] = 'None'
                mr_assigned[id].insert(0, 'None')

        # =============================================================================
        # FINAL DATA PROCESSING AND OUTPUT
        # =============================================================================
        
        # Filter out records with no valid breakpoint assignments
        mr_assigned2 = {}
        for id in mr_assigned:
            case = mr_assigned[id][2]  # Breakpoint classification
            if case != 'None':
                mr_assigned2[id] = mr_assigned[id]

        # Sort results by genomic coordinates
        sorted_dic_2 = dict(sorted(mr_assigned2.items(), key=lambda x: (x[1][0], x[1][1])))

        # Save final results to pickle
        try:
            with open(final_pickle_file, 'wb') as file:
                pickle.dump(sorted_dic_2, file)
            print(f"Pickle file created: {final_pickle_file}")
        except Exception as e:
            print(f"Error creating pickle file {final_pickle_file}: {e}")

        # Export final results to CSV
        try:
            with open(file_path, 'w', newline='') as csvfile:
                fieldnames = ['Key'] + [f'Value_{i+1}' for i in range(max(len(values) for values in sorted_dic_2.values()))]
                csv_writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
                csv_writer.writeheader()
                
                for key, values in sorted_dic_2.items():
                    row_data = {'Key': key}
                    for i, value in enumerate(values, start=1):
                        row_data[f'Value_{i}'] = value
                    csv_writer.writerow(row_data)
            print(f"CSV file created: {file_path}")
        except Exception as e:
            print(f"Error creating CSV file {file_path}: {e}")

    except Exception as e:
        print(f"Error processing file {picklepath}: {e}")

print("Genomic breakpoint analysis complete!")
