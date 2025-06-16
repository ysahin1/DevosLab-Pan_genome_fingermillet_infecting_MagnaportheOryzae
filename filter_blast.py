import os

# Define the input directory containing BLAST alignment result files
pathtoPot2LocInStrains = '../blastn_output_folder/'

# Define the output directory to save extracted read ID lists
output_path = '../idlist/'

# Collect all .txt files in the input directory
filelist = []
for file in os.listdir(pathtoPot2LocInStrains):
    if file.endswith('.txt'):
        filelist.append(os.path.join(pathtoPot2LocInStrains, file))

# Process each alignment file
for alignmentfilespath in filelist:
    # Extract filename
    alignmentfiles = alignmentfilespath.strip().split('/')[-1]

    # Parse metadata from filename assuming the pattern: param_outfmt_readname_strain.txt
    strain_name1 = alignmentfiles.strip().split('_')[3]
    strain_name = strain_name1.strip().split('.')[0]
    readname = alignmentfiles.strip().split('_')[2]
    outfmt = alignmentfiles.strip().split('_')[1]
    para = alignmentfiles.strip().split('_')[0]

    # Construct the output filename
    blastfile = f'{para}_{outfmt}_{readname}_idlist_{strain_name1}'

    # Open the output file for writing
    with open(f'{output_path}/{blastfile}', 'w') as out:
        # Open and read the alignment input file
        with open(alignmentfilespath) as file:
            for line in file:
                x = line.strip().split('\t')
                readid = x[1]  # Assuming second column contains read ID
                result = f'{readid}\n'
                out.write(result)

