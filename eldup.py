import os

# Function to filter rows based on breakpoint position and count status
def filter_table(data):
    filtered_data = []
    i = 0

    # Iterate through the list of data
    while i < len(data):
        # Check if the next row exists and has the same breakpoint (column 7)
        if i < len(data) - 1 and data[i][7] == data[i + 1][7]:
            # Both rows have "no_count" – skip one of them
            if data[i][13] == "no_count" and data[i + 1][13] == "no_count":
                i += 1
            # Both rows have "count" – keep both
            elif data[i][13] == "count" and data[i + 1][13] == "count":
                filtered_data.append(data[i])
                filtered_data.append(data[i + 1])
                i += 1
            # One "count", one "no_count" – keep the one with "count"
            elif data[i][13] != data[i + 1][13]:
                if data[i][13] == "count" and data[i + 1][13] == 'no_count':
                    filtered_data.append(data[i])
                elif data[i][13] == "no_count" and data[i + 1][13] == 'count':
                    filtered_data.append(data[i + 1])
                i += 2  # Move past both rows
        else:
            # If there's no pair to compare, just keep the current row
            filtered_data.append(data[i])
        i += 1

    return filtered_data

# Path to main directory with result files
pathToResTable = '../TotalBPinStrains'

reslist = []
# Collect all result files with .txt extension
for file in os.listdir(pathToResTable):
    if file.endswith('.txt'):
        reslist.append(os.path.join(pathToResTable, file))

# Process each result file
for respath in reslist:
    resfile = respath.strip().split('/')[-1]
    strain_name = resfile.strip().split('_')[1].split('.')[0]
    file_name = f'bpCount_{strain_name}.txt'

    data = []

    # Read original breakpoint data
    with open(respath, 'r') as file1:
        for line1 in file1:
            x1 = line1.strip().split('\t')
            data.append(x1)

    # Append breakpoint data from various TE families and directions
    sources = [
        "Grosshopper/start", "Grosshopper/end",
        "Fosbury/start", "Fosbury/end",
        "Pyret/start", "Pyret/end",
        "Family1/start", "Family1/end",
        "Family13/start", "Family13/end",
        "Inago2/start", "Inago2/end"
    ]

    for source in sources:
        filepath = f'/{source}/BpAssociated/TotalBPinStrains/{file_name}'
        with open(filepath, 'r') as f:
            for line in f:
                data.append(line.strip().split('\t'))
        print(f'{source} was added')

    # Sort data by column 12 (likely strain ID or chromosome) and column 7 (breakpoint position)
    data = sorted(data, key=lambda x: (x[12], x[7]))

    # Filter data to resolve breakpoint duplicates
    filtered_data = filter_table(data)

    # Additional labeling to mark breakpoints with no associated TE
    i = 0
    while i < len(filtered_data):
        if i < len(filtered_data) - 1 and filtered_data[i][0] == 'No_TE' or filtered_data[i][13] == 'count':
            i += 1
        elif i < len(filtered_data) - 1 and filtered_data[i][7] != filtered_data[i+1][7]:
            filtered_data[i][0] = 'No_TE'
            i += 1
        else:
            i += 1

    # Define output path
    output_filepath = f"../merged_bpcounts/merged_bpCounts_{strain_name}_.txt"

    # Save merged and filtered result
    with open(output_filepath, "w") as file:
        for row in filtered_data:
            file.write("\t".join(map(str, row)) + "\n")

    print(f"Data saved to {output_filepath}")

