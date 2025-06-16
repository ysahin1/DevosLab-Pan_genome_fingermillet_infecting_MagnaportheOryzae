import os
import pickle
import csv
import re
import logging
from Bio import SeqIO

# ---- Parameters and Paths ---- #
initial_parameter1 = 20
initial_parameter2 = 130
WantedStrain = []  # not used, could be removed

# Input/output paths
pathtobreakpoints = '../BP/breakpoints_locations.pkl'
distance = '../distance/distance.txt'
infopickle = 'pickle/'
csv_file_paths = 'infoCVS/'
pathtoBAMs = 'BAMs'
pathtoalignloc = 'LocInStrains'
Strainfasta = "Fastas/"

# ---- Load Data ---- #
# Load breakpoints
with open(pathtobreakpoints, 'rb') as file:
    breakpoints = pickle.load(file)

# Load distance data
with open(distance) as fff:
    distance_data = [line.strip().split('\t') for line in fff]

# Collect files
transloclist = [os.path.join(pathtoalignloc, f) for f in os.listdir(pathtoalignloc) if f.endswith('.txt')]
Strainfastalist = [os.path.join(Strainfasta, f) for f in os.listdir(Strainfasta) if f.endswith('.fasta')]
bamfilelist = [os.path.join(pathtoBAMs, f) for f in os.listdir(pathtoBAMs) if f.endswith('.bam')]

# ---- Helper Functions ---- #
def extract_numbers(part):
    return [int(num) for num in re.findall(r'\d+', part)]

def compare_numbers(s):
    parts = re.split(r'[MID]', s)
    left_numbers, right_numbers = [], []
    for i, part in enumerate(parts):
        nums = extract_numbers(part)
        (left_numbers if i == 0 else right_numbers).extend(nums)
    max_left = max(left_numbers) if left_numbers else float('-inf')
    max_right = max(right_numbers) if right_numbers else float('-inf')
    return [(s, 'left' if max_left > max_right else 'right', max(max_left, max_right))]

# ---- Main Processing Loop ---- #
trans_overlap_count_dic = {}
cigar_count_dic = {}
mapg_fail_dic = {}

for bam in bamfilelist:
    bamname = os.path.basename(bam).split('_')[2]
    readname = os.path.basename(bam).split('_')[1]
    strain_name = bamname.split('.')[0]
    entry_id = f'{strain_name}_{readname}'
    pathforpickle = f'{infopickle}info_{entry_id}.pkl'

    if os.path.exists(pathforpickle):
        logging.debug(f"File {pathforpickle} exists")
        continue

    try:
        print(f'Processing: {entry_id}')
        info = {}
        trans_overlap_count = 0
        cigar_count = 0
        mapg_fail = 0

        # Read BAM (text-based SAM-like format)
        with open(bam) as f:
            for line in f:
                x = line.strip().split('\t')
                id, chr, pos, mapg, cigar, seq = x[0], x[2], x[3], x[4], x[5], x[9]
                if int(mapg) < 0:
                    mapg_fail += 1
                if cigar == "*":
                    cigar_count += 1

                # Check if position overlaps breakpoints
                for loc, break_list in breakpoints.items():
                    for i in break_list:
                        strain = i[0]
                        if strain != strain_name:
                            continue
                        if len(i[1]) < 5:
                            continue

                        disup = int(i[1][3])
                        forward_pos = int(loc) - disup - 200
                        endSNP = i[1][0]
                        downward_pos = int(i[1][0])
                        downdis = int(i[1][4])
                        extend_down = downward_pos + downdis + 200

                        # Check overlap region (upstream)
                        if forward_pos <= int(pos) <= int(loc) and i[1][2] == chr:
                            trans_overlap_count += 1
                            info.setdefault(id, [pos, loc, chr, endSNP, cigar, 'None', seq])
                        # Check downstream extension
                        elif downward_pos <= int(pos) <= extend_down and i[1][2] == chr:
                            trans_overlap_count += 1
                            info.setdefault(id, [pos, loc, chr, downward_pos, cigar, 'None', seq])

        print(f'Closed BAM: {entry_id}')
        trans_overlap_count_dic[entry_id] = trans_overlap_count
        cigar_count_dic[entry_id] = cigar_count
        mapg_fail_dic[entry_id] = mapg_fail

        # Add distance to info
        for id in info:
            bpstart, bpchr = str(info[id][1]), str(info[id][2])
            for line in distance_data:
                if bpstart == line[1] and bpchr == line[0]:
                    info[id].insert(1, line[2])  # upstream distance
                    break
        print(f'Distance added: {entry_id}')

        for id in info:
            bpend, bpchr = str(info[id][4]), str(info[id][3])
            for i in range(len(distance_data)):
                if bpend == distance_data[i][1] and bpchr == distance_data[i][0]:
                    info[id].insert(3, distance_data[i+1][2])  # downstream distance
                    break

        print(f'Below distance added: {entry_id}')

        # Get matching sequence from FASTA
        usethisfasta = next((f for f in Strainfastalist if entry_id == '_'.join(os.path.basename(f).split('.')[:2])), None)
        sequence = {}
        with open(usethisfasta, 'r') as fasta:
            for record in SeqIO.parse(fasta, "fasta"):
                sequence[record.id.split()[0]] = str(record.seq)

        for id in info:
            if id in sequence:
                info[id].insert(0, sequence[id])

        print(f'Sequence added: {entry_id}')

        # Add alignment coordinates
        blastfilename = next((f for f in transloclist if entry_id in f), None)
        with open(blastfilename) as f:
            alignloc = [line.strip().split('\t') for line in f]

        for id in info:
            for line in alignloc:
                if id == line[1]:
                    info[id].insert(1, line[8])  # align start
                    info[id].insert(2, line[9])  # align end
                    break

        print(f'Align locations added: {entry_id}')

        # Analyze CIGAR string
        for id in info:
            cigar1, cigar2 = info[id][9], info[id][11]
            info[id].insert(10, compare_numbers(cigar1)[0][1] if cigar1 != '*' else 'None')
            info[id].insert(12, compare_numbers(cigar2)[0][1] if cigar2 != '*' else 'None')

        # Add TNE flag
        for id in info:
            tne = int(info[id][2])
            info[id].insert(3, 'Yes' if tne < initial_parameter1 or tne > initial_parameter2 else 'No')

        print(f'CIGAR and TNE added: {entry_id}')

        # ---- Save Pickle and CSV ---- #
        with open(pathforpickle, 'wb') as file:
            pickle.dump(info, file)
        print(f"Pickle saved: {pathforpickle}")

        csv_file_path = f'{csv_file_paths}bp_{entry_id}.csv'
        with open(csv_file_path, 'w', newline='') as csvfile:
            fieldnames = ['Key'] + [f'Value_{i+1}' for i in range(max(len(values) for values in info.values()))]
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            for key, values in info.items():
                row = {'Key': key}
                row.update({f'Value_{i+1}': v for i, v in enumerate(values)})
                writer.writerow(row)
        print(f"CSV saved: {csv_file_path}")

    except Exception as e:
        print(f"Error processing {entry_id}: {e}")

# ---- Summary ---- #
for strain in trans_overlap_count_dic:
    print(f'{trans_overlap_count_dic[strain]} overlaps found in {strain}')
for strain in cigar_count_dic:
    print(f'CIGAR failed in {cigar_count_dic[strain]} reads of {strain}')
for strain in mapg_fail_dic:
    print(f'MAPQ failed in {mapg_fail_dic[strain]} reads of {strain}')

