import os

# Path to directory containing read ID lists (output of previous step)
pathtoPot2LocInStrainsIdlist = 'idlist/'

# Path to directory containing BAM files for different strains
bam_path = 'Path_to_BAM'

# Output directory where filtered BAM entries will be written
output_bam_filepath = '../BAMs'

# Collect all ID list files (assumed to be .txt files)
idfilelist = []
for file in os.listdir(pathtoPot2LocInStrainsIdlist):
    if file.endswith('.txt'):
        idfilelist.append(os.path.join(pathtoPot2LocInStrainsIdlist, file))

# Collect all BAM files from input BAM path
bamfilelist = []
for file in os.listdir(bam_path):
    if file.endswith('.bam'):
        bamfilelist.append(os.path.join(bam_path, file))

# For each ID list file, extract strain and construct output name
for idlistpath in idfilelist:
    # Extract metadata from filename
    alignmentfiles = idlistpath.strip().split('/')[-1]
    strain_name1 = alignmentfiles.strip().split('_')[4]
    strain_name = strain_name1.strip().split('.')[0]
    readname = alignmentfiles.strip().split('_')[2]
    outfmt = alignmentfiles.strip().split('_')[1]
    para = alignmentfiles.strip().split('_')[0]

    # Construct the desired output BAM filename
    bamfile_name = f'{para}_{readname}_{strain_name}.bam'

    # Find the matching BAM file for the current strain
    for bamfile in bamfilelist:
        bamfilesplit = bamfile.strip().split('/')[-1]
        bam_strain_name = bamfilesplit.strip().split('.')[0]

        # If strain names match, extract matching reads using samtools and fgrep
        if bam_strain_name == strain_name:
            os.system(
                f'samtools view {bamfile} | fgrep -f {idlistpath} > {output_bam_filepath}/{bamfile_name}'
            )

