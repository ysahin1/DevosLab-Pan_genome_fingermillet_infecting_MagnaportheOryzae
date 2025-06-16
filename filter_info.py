import pickle
import csv
import os

# === Step 1: Define paths ===
# Path to the folder containing the input .pkl files
pathtopicklefolder = '../pickle/'

# List to store full paths to .pkl files
picklelist = []
for file in os.listdir(pathtopicklefolder):
    if file.endswith('.pkl'):
        picklelist.append(os.path.join(pathtopicklefolder, file))

# Path to store filtered output .pkl files
pathtoOutputPickle = '../filteredpickle/'

# Path to store output CSV files
csv_file_path = '../info2CSV/'

# === Step 2: Process each .pkl file ===
for picklepath in picklelist:
    # Extract filename parts to determine strain and read name
    picklefile = picklepath.strip().split('/')[-1]
    readname1 = picklefile.strip().split('_')[2]
    readname = readname1.strip().split('.')[0]
    strain_name = picklefile.strip().split('_')[1]
    picklefile = f'{strain_name}_{readname}'

    # Load the .pkl file into a dictionary called 'info'
    with open(picklepath, 'rb') as file:
        info = pickle.load(file)

    # === Step 3: Filter records ===
    # Create a new dictionary to store filtered results
    info2 = {}
    for id in info:
        condition = info[id][3]     # Assumes this index indicates breakpoint proximity
        cigar1 = info[id][10]       # First CIGAR string
        cigar2 = info[id][12]       # Second CIGAR string

        # Keep only records marked 'Yes' and containing insertions (I) or deletions (D) in either CIGAR
        if condition == 'Yes':
            if "I" in cigar1 or "D" in cigar1:
                info2[id] = info[id]
            else:
                if "I" in cigar2 or "D" in cigar2:
                    info2[id] = info[id]

    # === Step 4: Save filtered dictionary to new .pkl file ===
    try:
        with open(f'{pathtoOutputPickle}/info2_{picklefile}.pkl', 'wb') as file:
            pickle.dump(info2, file)
    except Exception as e:
        print(f"Error creating pickle file {picklefile}: {e}")

    # === Step 5: Export filtered records to CSV ===
    try:
        with open(f'{csv_file_path}/bp2_{picklefile}.csv', 'w', newline='') as csvfile:
            # Prepare header with 'Key' and Value_1, Value_2, ..., based on longest value list
            fieldnames = ['Key'] + [f'Value_{i+1}' for i in range(max(len(values) for values in info2.values()))]

            csv_writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

            # Write the header row
            csv_writer.writeheader()

            # Write each entry in the dictionary
            for key, values in info2.items():
                row_data = {'Key': key}
                for i, value in enumerate(values, start=1):
                    row_data[f'Value_{i}'] = value
                csv_writer.writerow(row_data)
    except Exception as e:
        print(f"Error creating CSV file {csv_file_path}/bp2_{picklefile}.csv: {e}")

