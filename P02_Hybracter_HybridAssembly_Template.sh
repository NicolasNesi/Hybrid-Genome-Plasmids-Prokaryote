#!/bin/bash
#SBATCH -J P02_Hybracter                    # Job name 
#SBATCH -o P02_Hybracter.%a.%A.out          # File to which stdout will be written
#SBATCH -e P02_Hybracter.%a.%A.err          # File to which stderr will be written
#SBATCH -p court                            # Partition
#SBATCH -N 1                                # Number of node
#SBATCH -c 96                               # Number of cores/cpus
#SBATCH -t 02-00:00                         # Runtime in DD-HH:MM
#SBATCH --mem 400G                          # Memory for all cores in Gbytes
#SBATCH --mail-type=ALL                     # BEGIN,END,ALL
#SBATCH --mail-user=mymail                  # Email address

# ---------------------------------
# Version v1.9
# ---------------------------------
# authors: Nicolas Nesi
# University Caen Normandy, DYNAMICURE INSERM UMR 1311
# Date:27/06/2025
# ---------------------------------

# ---------------------------------
# To install a singurality image of hybracter run:
# srun --partition=2tcourt --mem=50G singularity pull --dir $IMAGE_DIR docker://quay.io/gbouras13/hybracter:x.x.x
# ---------------------------------

# ---------------------------------
# Variables and Paths
# ---------------------------------
# Positional parameters
GROUP="${1}"
PROJECT="${2}"
RUN="${3}"
MINLENGTH="${4}"
QSCOREONT="${5}"
MEDAKA="${6}"
LOGIC="${7}"
MINDEPTH="${8}"
SUBSAMPLE="${9}"
FILTER=5000

# Paths
FOLDER="/dlocal/home/2019013/Data/$GROUP/$PROJECT"
PATH1="$FOLDER/RawReads_Nanopore/$RUN"
PATH2="$FOLDER/RawReads_Illumina/$RUN"
PATH3="$FOLDER/Hybracter/$RUN"
IMAGE_DIR="/soft/2019013/Logiciels/hybracter"

# Image Singularity
containerImage="$IMAGE_DIR/hybracter_0.11.2.sif"

# Databases
BAKTADB="/dlocal/home/2019013/Databases/Bakta/db"
AMR="/dlocal/home/2019013/Databases/AMRFinderPlus/2024-01-31.1/"

RUNDATE="$(date '+%Y-%m-%d')"

# Threads
export THREADS="${SLURM_CPUS_PER_TASK}"
echo "number of threads used ${THREADS}"
# ---------------------------------

# ---------------------------------
# Make subdirectories
# ---------------------------------
for dir in \
	"${FOLDER}/CleanedReads_Nanopore/${RUN}" \
	"${FOLDER}/CleanedReads_Illumina/${RUN}" \
	"${FOLDER}/Hybracter/${RUN}" \
    "${FOLDER}/Plasmids/${RUN}" \
    "${FOLDER}/Chromosome/${RUN}" \
    "${FOLDER}/Bakta/${RUN}" \
    "${FOLDER}/Staramr/${RUN}" \
    "${FOLDER}/QC_Stats/${RUN}/PycoQC" \
    "${FOLDER}/QC_Stats/${RUN}/QUAST" \
	"${FOLDER}/QC_Stats/${RUN}/BUSCO" \
	"${FOLDER}/QC_Stats/${RUN}/Kraken2" \
	"${FOLDER}/QC_Stats/${RUN}/Seqkit" \
	"${FOLDER}/QC_Stats/${RUN}/gfa" \
	"${FOLDER}/QC_Stats/${RUN}/RawReads_Nanopore" \
	"${FOLDER}/QC_Stats/${RUN}/RawReads_Illumina" \
	"${FOLDER}/QC_Stats/${RUN}/CleanedReads_Nanopore" \
	"${FOLDER}/QC_Stats/${RUN}/CleanedReads_Illumina" \
    "${FOLDER}/Log/${RUN}_run${RUNDATE}" \
    "${FOLDER}/Scripts/Sbatch_log"
do
    if [ ! -d "$dir" ]; then
        mkdir -p "$dir"
        echo "Created directory: $dir"
    else
        echo "Directory already exists: $dir"
    fi
done
# ---------------------------------

# ---------------------------------
# Copy input data to LOCAL_WORK_DIR
# ---------------------------------
cp "${PATH1}"/*.fastq.gz "${LOCAL_WORK_DIR}"
cp "${PATH2}"/*_R*.fastq.gz "${LOCAL_WORK_DIR}"
cp "${FOLDER}/Scripts/hybracter_${RUN}.csv" "${LOCAL_WORK_DIR}"
cd "${LOCAL_WORK_DIR}"/

echo Analysis run the: "${RUNDATE}"
echo Working directory: "${PWD}"
echo Submitted with Variables: Group "${GROUP}" Project "${PROJECT}" Run "${RUN}" Minimum length Nanopore reads "${MINLENGTH}" Minimum Qscore Nanopore reads "${QSCOREONT}"
echo                           Medaka model "${MEDAKA}" Polishing selection results "${LOGIC}" Minimum depth "${MINDEPTH}" Subsampling depth "${SUBSAMPLE}"
# ---------------------------------

# ---------------------------------
# Hybracter run
# ---------------------------------
echo "Start Hybracter"

singularity exec "$containerImage" \
hybracter hybrid \
--input hybracter_"${RUN}".csv \
--datadir "${LOCAL_WORK_DIR}" \
--contaminants lambda \
--min_length "${MINLENGTH}" \
--min_quality "${QSCOREONT}" \
--subsample_depth "${SUBSAMPLE}" \
--min_depth "${MINDEPTH}" \
--depth_filter 0.25 \
--flyeModel "--nano-hq" \
--logic "${LOGIC}" \
--no_medaka \
--threads "${THREADS}" \
--output "${RUN}"_output/

# Move hybracter outputs
cp -R "${RUN}"_output/* "${PATH3}"
cp "${PATH3}"/hybracter.log "${FOLDER}/Log/${RUN}_run${RUNDATE}/hybracter_run${RUNDATE}.log"
# ---------------------------------

# ---------------------------------
# Assembly Annotation
# ---------------------------------
module purge
module load py_env/anaconda3
eval "$(conda shell.bash hook)"
conda activate bakta_env
# include:
# bakta=1.11.3

echo "Annotation of complete assemblies"

for file in "${PATH3}"/FINAL_OUTPUT/complete/*_final.fasta; do
    base=$(basename "${file}" _final.fasta)

    bakta \
        --db "${BAKTADB}" \
        "${file}" \
        --prefix "${base}" \
        --threads "${THREADS}" \
		--force \
		--verbose \
        --output "${FOLDER}/Bakta/${RUN}/Complete/${base}" \
        2>&1 | tee -a "${FOLDER}/Log/${RUN}_run${RUNDATE}/${base}_bakta_hybracter_run${RUNDATE}.log"
done
# ---------------------------------

# ---------------------------------
# Plasmid Typing
# ---------------------------------
conda deactivate
conda activate staramr_env
# include:
# staramr=0.11.0

echo "Perform ResFinder and PlasmidFinder on full assemblies"

staramr search \
--nprocs "${THREADS}" \
--output-dir out/ \
"${PATH3}"/FINAL_OUTPUT/complete/*_final.fasta 2>&1 | tee -a "${FOLDER}/Log/${RUN}_run${RUNDATE}/staramr_run${RUNDATE}.log"

mkdir -p "${FOLDER}/Staramr/${RUN}"
mv out/* "${FOLDER}/Staramr/${RUN}"
# ---------------------------------

# ---------------------------------
# Antimicrobial Resistance Gene Finder
# ---------------------------------
conda deactivate
conda activate bakta_env
# include:
# amrfinderplus=4.2.7

echo "Research of AMR genes in plasmids"

for file in "${PATH3}"/FINAL_OUTPUT/complete/*_plasmid.fasta; do
    base=$(basename "${file}" _plasmid.fasta)
    amrfinder \
    --database "${AMR}" \
    --nucleotide "${file}" \
    --output "${FOLDER}/Plasmids/${RUN}/${base}_AMRFinder_Plasmids_Hybracter_report.txt" \
    --nucleotide_output "${FOLDER}/Plasmids/${RUN}/${base}_AMRfinder_Plasmids_Hybracter.fasta" \
    --plus \
    --threads "${THREADS}" \
    --name "${base}" 2>&1 | tee -a "${FOLDER}/Log/${RUN}_run${RUNDATE}/${base}_AMRfinder_Plasmids_Hybracter_run${RUNDATE}.log"
done
# ---------------------------------

# ---------------------------------
# OXA-48 Results
# ---------------------------------
# Check if report files exist
shopt -s nullglob
reports=("${FOLDER}/Plasmids/${RUN}"/*_AMRFinder_Plasmids_Hybracter_report.txt)
if [ ${#reports[@]} -eq 0 ]; then
    echo "No AMRFinder report files found in ${FOLDER}/Plasmids/${RUN}."
    exit 1
fi

# Extract lines containing 'OXA-48'
cd "${FOLDER}/Plasmids/${RUN}"
for report in "${reports[@]}"; do
    basename_report=$(basename "${report}")
    grep "OXA-48" "${report}" | sed "s|^|${basename_report}: |" >> "List_OXA-48_Plasmids_Hybracter_${RUN}.txt"
done
# ---------------------------------

# ---------------------------------
# Copy Hybracter outputs
# ---------------------------------

conda deactivate
conda activate NN_env
# include:
# seqkit=2.9.0

echo "Copy statistic tsv files of Hybracter and plassembler"

mkdir -p "${FOLDER}/QC_Stats/${RUN}/Assemblies/"
cp "${PATH3}/FINAL_OUTPUT/hybracter_summary.tsv" "${FOLDER}/QC_Stats/${RUN}/Assemblies/"
cp "${PATH3}/supplementary_results/plassembler_all_assembly_summary/plassembler_assembly_info.tsv" "${FOLDER}/QC_Stats/${RUN}/Assemblies/"

echo "Copy plasmids contigs"

cp "${PATH3}"/FINAL_OUTPUT/complete/*_plasmid.fasta "${FOLDER}/Plasmids/${RUN}/"

for f in "$FOLDER/Plasmids/${RUN}"/*.fasta; do 
    seqkit seq -w0 "$f" > "${f}.tmp" && mv "${f}.tmp" "$f"
done

for f in "$FOLDER/Plasmids/${RUN}"/*_plasmid.fasta; do 
    grep -A 1 "circular=true" "$f" > "${f%.fasta}"_circular.fasta
done

find "$FOLDER/Plasmids/${RUN}"/ -type f -size 0 -delete

echo "Copy chromosome assemblies"

cp "${PATH3}"/FINAL_OUTPUT/complete/*_chromosome.fasta "${FOLDER}/Chromosome/${RUN}/"

for f in "${FOLDER}/Chromosome/${RUN}"/*.fasta; do 
    seqkit seq -w0 "$f" > "${f}.tmp" && mv "${f}.tmp" "$f"
done

echo "Copy cleaned Illumina and Nanopore reads and Fastp et Seqkit reports"

cp "${PATH3}"/processing/qc/fastp/*.fastq.gz "${FOLDER}/CleanedReads_Illumina/${RUN}/"
rename "_1" "_cleaned_R1" "$FOLDER/CleanedReads_Illumina/$RUN"/*_1.fastq.gz
rename "_2" "_cleaned_R2" "$FOLDER/CleanedReads_Illumina/$RUN"/*_2.fastq.gz

cp "${PATH3}"/processing/qc/fastp/*.html "${FOLDER}/QC_Stats/${RUN}/CleanedReads_Illumina/"
rename ".html" "_Illumina_ReportFastp.html" "${FOLDER}/QC_Stats/${RUN}/CleanedReads_Illumina"/*.html
cp "${FOLDER}/QC_Stats/${RUN}/CleanedReads_Illumina"/*.html "${FOLDER}/QC_Stats/${RUN}/RawReads_Illumina/"

cp "${PATH3}"/processing/qc/*_filt_trim.fastq.gz "${FOLDER}/CleanedReads_Nanopore/${RUN}/"
rename "_filt_trim" "_cleaned_min$MINLENGTH" "${FOLDER}/CleanedReads_Nanopore/${RUN}"/*.fastq.gz

cp "${PATH3}"/processing/qc/seqkit/*.txt "${FOLDER}/QC_Stats/${RUN}/Seqkit/"
rename "_long" "_Nanopore" "${FOLDER}/QC_Stats/${RUN}/Seqkit"/*_long.txt
rename "_r1" "_Illumina_R1" "${FOLDER}/QC_Stats/${RUN}/Seqkit"/*_r1.txt
rename "_r2" "_Illumina_R2" "${FOLDER}/QC_Stats/${RUN}/Seqkit"/*_r2.txt
sed -i "s|${RUN}_output/processing/qc/fastp/||g" "${FOLDER}/QC_Stats/${RUN}/Seqkit/"*.txt
# ---------------------------------

# ---------------------------------
# Filter Plasmids contigs and splits
# ---------------------------------
echo "Filter plasmids contigs by length and split fasta"

for f in "$FOLDER/Plasmids/${RUN}"/*_plasmid.fasta; do
  seqkit seq -m "${FILTER}" "$f" > "${f%.fasta}"_contigssup"${FILTER}".fasta
done

for f in "$FOLDER/Plasmids/${RUN}"/*_contigssup"${FILTER}".fasta; do
  seqkit split -i -O SplitFasta "$f"
done

for f in "$FOLDER/Plasmids/${RUN}"/SplitFasta/*_plasmid_contigssup"$FILTER".part_*.fasta; do
      rename "_plasmid_contigssup5000.part" "" "$f"
done
# ---------------------------------

# ---------------------------------
# Each plasmids Annotation
# ---------------------------------
conda deactivate
conda activate bakta_env
# include:
# bakta=1.11.3

echo "Annotation of circularized plasmids for sample: $SAMPLE"

for file in "$FOLDER/Plasmids/${RUN}"/SplitFasta/*.fasta; do
    base=$(basename "${file}" .fasta)

    bakta \
        --db "${BAKTADB}" \
        "${file}" \
        --prefix "${base}" \
        --threads "${THREADS}" \
		--force \
		--verbose \
        --output "${FOLDER}/Bakta/${RUN}/Plasmids/${base}" \
        2>&1 | tee -a "${FOLDER}/Log/${RUN}_run${RUNDATE}/${base}_bakta_plasmids_run${RUNDATE}.log"
done
# ---------------------------------

# ---------------------------------
# Cleaning Working Directory
# ---------------------------------
rm -R "${LOCAL_WORK_DIR}"/*
# ---------------------------------

# ---------------------------------
mv "${FOLDER}"/Scripts/P02_Hybracter.* "${FOLDER}"/Scripts/Sbatch_log
# ---------------------------------
