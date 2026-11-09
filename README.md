# Transposon-linked deletions in *Magnaporthe oryzae* strains that infect finger millet

Scripts (Devos Lab) that find structural breakpoints (deletions) in *Magnaporthe oryzae* strains
and check whether a transposable element (TE) sits at the breakpoint. Short reads from each strain are
searched for a TE sequence with BLAST. Reads that contain the TE and map close to a known breakpoint
are then classified and counted per strain and per chromosome.

## What you need

- Python 3 with [Biopython](https://biopython.org/)
- [BLAST+](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html) (`makeblastdb`, `blastn`)
- [samtools](http://www.htslib.org/) and [seqtk](https://github.com/lh3/seqtk)

Inputs that the scripts expect (the paths are set at the top of each script):

| Input | Description |
|---|---|
| Trimmed reads (`*_trimed.fastq`) | Paired-end reads for each strain |
| BAM files | Reads aligned to the reference genome, one file per strain |
| `transposon.fa` | Sequence of the TE to search for (for example Pot2) |
| `breakpoints_locations.pkl` | Known breakpoints: a dictionary keyed by breakpoint start; each value lists the strain, end position, size and chromosome |
| `distance.txt` | Tab-separated: chromosome, SNP position, distance to the next SNP |

## Pipeline, in run order

| Step | Script | What it does | Writes to |
|---|---|---|---|
| 1 | `cpfasta.sh` | Converts each trimmed FASTQ to FASTA with `seqtk` and renames it `<strain>.R<1/2>.fasta` | FASTA folder |
| 2 | `mkdbblast.sh` | Builds a BLAST database from each strain's reads and searches them with the TE sequence (`blastn-short`, tabular output) | `Pot2LocInStrains/` |
| 3 | `filter_blast.py` | Takes the read IDs of the BLAST hits | `idlist/` |
| 4 | `filter_bam.py` | Pulls those reads out of each strain's BAM (`samtools view \| fgrep -f`) | `BAMs/` (SAM text) |
| 5 | `bp.py` | Keeps reads that map within a breakpoint's flanking window (SNP distance + 200 bp). Adds the read sequence, the TE alignment coordinates, the flanking SNP distances and the CIGAR side with the longest match. Flags reads whose TE hit starts before base 20 or after base 130 (TE at the read end) | `pickle/`, `infoCVS/` |
| 6 | `filter_info.py` | Keeps flagged reads that also have an insertion or deletion (`I`/`D`) in their CIGAR | `filteredpickle/`, `info2CSV/` |
| 7 | `CaseAssigner.py` | Merges R1 and R2 for each read pair and assigns a breakpoint "case" by read orientation and whether the mate maps | `mergedpickle/`, `infoResult/` |
| 8 | `CountTans.py` | Groups consecutive reads into intervals and counts, per strain and chromosome (Chr01–Chr07), whether a TE is supported at one end or at both ends | `CountsInEachStrain/all_strainAssocount_result.csv` |
| 9 | `CountBPs.py` | Writes one table per strain: each breakpoint with its supporting read, case, flanking SNP distances and SNP count. A breakpoint is marked `count` when the flank is shorter than 3,000 bp. Duplicate rows are collapsed with a count | `TotalBPinStrains/bpCount_<strain>.txt` |
| 10 | `eldup.py` | Merges the per-strain tables from all TE families (Grasshopper, Fosbury, Pyret, Family1, Family13, Inago2; start and end) and removes duplicate breakpoints, keeping `count` rows over `no_count` rows | per-strain table |

## Notes

- The paths are hard-coded relative folders (for example `../BP/`, `../distance/`). Edit them before running.
  Some later steps read folders with a `2` suffix (`mergedpickle2/`, `infoResult2/`), so rename or copy the
  earlier outputs to match.
- Strain and read names come from the file names (the parts split on `_` and `.`). Keep the naming pattern
  from steps 1–2.
- Steps 2–10 are run once for each TE family; `eldup.py` then combines the families.
