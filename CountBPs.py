import pickle
import os

# File paths
pathTobpdic = '../BP/breakpoints_locations.pkl'  # Input pickle file with breakpoint data
pathToResTable = '../infoResult2'  # Directory containing CSV files for strains
output_path = '../TotalBPinStrains'  # Output folder for results
distance = '../distance/distance.txt'  # File containing SNP distances
transposon_type = ''  # Label used for all output lines

# Read distance file into memory
with open(distance) as fff:
    distance_data = [line.strip('\n').strip().split('\t') for line in fff]

# Load breakpoint dictionary
with open(pathTobpdic, 'rb') as file:
    bpdic = pickle.load(file)

countdic = {}  # Stores breakpoints grouped by strain

# Populate countdic with data from bpdic
for bp in bpdic:
    for i in bpdic[bp]:
        strain = i[0]
        bpend = i[1][0]
        bpsize = i[1][1]
        bpchr = i[1][2]
        if strain not in countdic:
            countdic[strain] = []
        countdic[strain].append([bp, bpend, bpsize, bpchr])

# Match each bp start to the SNP distance data
for strain in countdic.keys():
    for idx, i in enumerate(countdic[strain]):
        bpchr = i[3]
        bpstart = i[0]
        for line_idx, line in enumerate(distance_data):
            dischr = str(line[0])
            bp1 = line[1]
            dis = line[2]
            if int(bpstart) == int(bp1) and str(bpchr) == str(dischr):
                countdic[strain][idx].insert(3, dis)        # Insert left SNP distance
                countdic[strain][idx].append(line_idx + 1)  # Add index of match
                break

# Match each bp end to SNP distance data
for strain in countdic.keys():
    for idx, i in enumerate(countdic[strain]):
        bpchr = i[4]
        bpend = i[1]
        for ii in range(len(distance_data)):
            dischr = str(distance_data[ii][0])
            bp1 = distance_data[ii][1]
            dis = distance_data[ii][2]
            if int(bpend) == int(bp1) and bpchr == dischr:
                if ii + 1 < len(distance_data):
                    countdic[strain][idx].insert(4, distance_data[ii + 1][2])  # Insert right SNP distance
                    countdic[strain][idx].append(ii + 1)  # Add right index
                    break

reslist = []  # List of all result table file paths

# Collect all result CSVs
for file in os.listdir(pathToResTable):
    if file.endswith('.csv'):
        reslist.append(os.path.join(pathToResTable, file))

# Process each result table file
for respath in reslist:
    resfile = respath.strip().split('/')[-1]
    strain_name = resfile.strip().split('_')[2]
    file_name = f'bpCount_{strain_name}.txt'
    strainWanted = strain_name

    # Read CSV content into info dictionary
    info = {}
    with open(respath) as file2:
        for line in file2:
            if 'Key' in line:
                continue
            x = line.strip().split(',')
            id = x[0]
            if id not in info:
                info[id] = []
            info[id].append(x[2])   # actual location
            info[id].append(x[3])   # case
            info[id].append(x[9])   # column 9
            info[id].append(x[10])  # column 10
            info[id].append(x[11])  # start
            info[id].append(x[12])  # unused?
            info[id].append(x[14])  # end

    # Write processed data to text file
    with open(f'{output_path}/{file_name}', 'w') as out:
        for strain in countdic:
            if strain == strainWanted:
                for pos in countdic[strain]:
                    found = False
                    for id in info:
                        bpstart = info[id][4]
                        bpend = info[id][6]
                        if int(pos[0]) == int(bpstart):
                            snpcount = int(pos[7]) - int(pos[6]) + 1
                            calvalue = (int(pos[1]) + int(pos[4])) - (int(pos[0]) - int(pos[3]))

                            # Logic for overlap and size threshold
                            if int(info[id][0]) < int(pos[0]):
                                result = (f"{transposon_type}\t{id}\t{info[id][1]}\t{calvalue}\t{info[id][0]}\t{info[id][2]}"
                                          f"\t{pos[3]}\t{pos[0]}\t{pos[4]}\t{pos[1]}\t{pos[2]}\t{snpcount}\t{pos[5]}")
                                result += "\tcount\n" if int(pos[3]) < 3000 else "\tno_count\n"
                            elif int(info[id][0]) > int(pos[1]):
                                result = (f"{transposon_type}\t{id}\t{info[id][1]}\t{calvalue}\t{info[id][0]}\t{info[id][2]}"
                                          f"\t{pos[3]}\t{pos[0]}\t{pos[4]}\t{pos[1]}\t{pos[2]}\t{snpcount}\t{pos[5]}")
                                result += "\tcount\n" if int(pos[4]) < 3000 else "\tno_count\n"
                            out.write(result)
                            found = True
                    # If no matching ID found, write partial data
                    if not found and len(pos) > 7:
                        snpcount = int(pos[7]) - int(pos[6]) + 1
                        result = f"{transposon_type}\t\t\t\t\t\t{pos[3]}\t{pos[0]}\t{pos[4]}\t{pos[1]}\t{pos[2]}\t{snpcount}\t{pos[5]}\tno_count\n"
                        out.write(result)

    print(f"TXT file created: {output_path}/{file_name}")

    # Post-process file to count duplicate entries by position
    data = []
    with open(f'{output_path}/{file_name}', 'r') as file1:
        for line1 in file1:
            x1 = line1.strip().split('\t')
            data.append(x1)

    column_index = 7  # Position used for duplicate tracking

    # Dictionary to hold and count unique rows
    unique_data = {}
    duplicate_counts = {}

    for row in data:
        key = row[column_index]
        if key in unique_data:
            duplicate_counts[key] = duplicate_counts.get(key, 1) + 1
        else:
            unique_data[key] = row
            duplicate_counts[key] = 1

    # Add duplicate count to each unique row
    result = []
    for key, row in unique_data.items():
        row_with_count = row + [duplicate_counts[key]]
        result.append(row_with_count)

    # Save the deduplicated result with counts
    output_filepath = f'{output_path}/{file_name}'
    with open(output_filepath, "w") as file:
        for row in result:
            file.write("\t".join(map(str, row)) + "\n")

    print(f"Data saved to {output_filepath}")

