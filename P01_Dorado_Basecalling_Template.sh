#!/bin/bash
#SBATCH -J P01_DoradoGPU                          # Job name
#SBATCH -o P01_DoradoGPU.%a.%A.out                # File to which stdout will be written
#SBATCH -e P01_DoradoGPU.%a.%A.err                # File to which stderr will be written
#SBATCH -p gpu_all                                # Partition
#SBATCH -N 1                                      # Number of node - 1 Milan node has 32-core, 500 Go RAM and GPU Ampere A100
#SBATCH --gres gpu:8                              # Number GPU
#SBATCH --ntasks-per-node 1                       # Number of cores/cpus
#SBATCH -c 28                                     # Number of cpus per task
#SBATCH -t 00-08:00                               # Runtime in DD-HH:MM
#SBATCH --mem 100G                                # Memory for all cores in Gbytes
#SBATCH --mail-type=ALL                           # BEGIN,END,ALL
#SBATCH --mail-user=mymail                        # Email address

# ---------------------------------
# Version v1.9
# ---------------------------------
# author: Nicolas Nesi
# University Caen Normandy, DYNAMICURE INSERM UMR 1311
# Date: 27/06/2024
# ---------------------------------

# ---------------------------------
# Variables and Paths
# ---------------------------------
# Positional parameters
GROUP="${1}"
PROJECT="${2}"
RUN="${3}"
MODEL="${4}"
BARCODE="${5}"
QSCOREDORADO="${6}"

# Paths
FOLDER="/dlocal/home/2019013/Data/$GROUP/$PROJECT"
SCRIPTS="$FOLDER/Scripts"
PATH2="$FOLDER/RawReads_Nanopore/$RUN"

RUNDATE="$(date '+%Y-%m-%d')"

# Threads
export THREADS="${SLURM_CPUS_PER_TASK}"
echo "number of threads used ${THREADS}"
# ---------------------------------

# ---------------------------------
# Copy input data to LOCAL_WORK_DIR
# ---------------------------------
cp -R "${FOLDER}"/Nanopore_outputRun/"${RUN}"/no_sample_id/*/pod5_skip "${LOCAL_WORK_DIR}"
cd "${LOCAL_WORK_DIR}"/
mkdir demultiplex/

echo Analysis run the: "${RUNDATE}"
echo Working directory: "${PWD}"
echo Submitted with variables: Group "${GROUP}" Project "${PROJECT}" Run "${RUN}" Dorado-model "${MODEL}" NanoporeKit "${BARCODE}" Qscore "${QSCOREDORADO}"
# ---------------------------------

# ---------------------------------
# Make subdirectories
# ---------------------------------
for dir in \
    "${FOLDER}/Dorado_Basecalling/${RUN}" \
    "${FOLDER}/RawReads_Nanopore/${RUN}" \
    "${FOLDER}/RawReads_Illumina/${RUN}" \
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
# Modules
# ---------------------------------
module purge
module load bio/dorado/0.9.6
# ---------------------------------

# ---------------------------------
# Dorado GPU
# ---------------------------------
echo "Running Basecalling, Demultiplexing and Adapter/barcode trimming using Dorado"

dorado basecaller \
/soft/2019013/Logiciels/Dorado-models/"${MODEL}" \
"${LOCAL_WORK_DIR}"/pod5_skip/ \
--device "cuda:all" \
--min-qscore "${QSCOREDORADO}" \
--chunksize 10000 \
--no-trim \
--verbose > "${RUN}".basecalls.bam

dorado summary "${RUN}".basecalls.bam > demultiplex/sequencing_summary.txt

dorado demux \
--kit-name "${BARCODE}" \
--threads "${THREADS}" \
--sample-sheet "${SCRIPTS}"/Sample_sheet_"${RUN}".csv \
--emit-summary \
--output-dir demultiplex/ "${RUN}".basecalls.bam

for f in demultiplex/*.bam; do
   base=$(basename "$f" .bam)
   dorado trim \
   --no-trim-primers \
   --sequencing-kit "${BARCODE}" \
   --threads "${THREADS}" \
   "$f" > demultiplex/"${base}"-trimmed.bam
done
# ---------------------------------

# ---------------------------------
# Environments
# ---------------------------------
module purge
module load py_env/anaconda3
eval "$(conda shell.bash hook)"
conda activate clean-ONT_env
# include:
# samtools=1.18
# ---------------------------------

# ---------------------------------
# Convert uBAM to Fastq
# ---------------------------------
echo "Convert uBAM to Fastq"

for f in demultiplex/*-trimmed.bam; do
 base=$(basename "$f" -trimmed.bam)
 samtools bam2fq -@ "${THREADS}" "$f" > demultiplex/"${base}".fastq
done
# ---------------------------------

# ---------------------------------
# Move and rename fastq
# ---------------------------------
echo "Move and rename fastq"

mv demultiplex/ "${FOLDER}"/Dorado_Basecalling/"${RUN}"

mv "${RUN}".basecalls.bam "${FOLDER}"/Dorado_Basecalling/"${RUN}"

cd "${FOLDER}"/Dorado_Basecalling/"${RUN}"/demultiplex/

mv *_unclassified.fastq NoBarcode.fastq

for f in *.fastq; do
  new_name=$(echo "$f" | sed 's/.*_//')
  mv "$f" "$new_name"
done

mv "${FOLDER}"/Dorado_Basecalling/"${RUN}"/demultiplex/*.fastq "${PATH2}"/

pigz -f "${PATH2}"/*.fastq
# ---------------------------------

# ---------------------------------
# Copy Illumina Raw Reads
# ---------------------------------
echo "Copying Illumina raw reads from Illumina_outputRun/ to RawReads_Illumina/"

# Check if the target directory already contains .fastq.gz files
if ! compgen -G "${FOLDER}/RawReads_Illumina/${RUN}/*.fastq.gz" > /dev/null; then
    # Attempt to copy files
    if cp "${FOLDER}"/Illumina_outputRun/"${RUN}"/Alignment_1/*/Fastq/*.fastq.gz "${FOLDER}/RawReads_Illumina/${RUN}/"; then
        echo "Copy completed successfully."
    else
        echo "Error: Copy operation failed. Attempting to move files instead."

        # Loop through the list of sample names and move files
        while read -r line; do
            mv "${FOLDER}/RawReads_Illumina/${line}"_R*.fastq.gz "${FOLDER}/RawReads_Illumina/${RUN}/"
        done < "${FOLDER}/Scripts/List_samples_${RUN}.txt"

        echo "Move operation completed."
    fi
else
    echo "Fastq files already exist in the target directory."
fi
# ---------------------------------

# ---------------------------------
# Environments
# ---------------------------------
conda deactivate
conda activate QC-ONT_env
# Include:
# PycoQC=2.5.2
# ---------------------------------

# ---------------------------------
# Quality Control with PycoQC
# ---------------------------------
echo "QC with PycoQC"

cd "${FOLDER}"/Dorado_Basecalling/"${RUN}"/demultiplex/

# format summary files for pycoQC
{ head -n 1 sequencing_summary.txt; tail -n +2 sequencing_summary.txt | sort -k2,2; } > sequencing_summary_sorted.txt
{ head -n 1 barcoding_summary.txt; tail -n +2 barcoding_summary.txt | sort -k2,2; } > barcoding_summary_sorted.txt
header=$(head -n 1 sequencing_summary_sorted.txt)
awk 'NR==FNR {barcodes[FNR]=$NF; next} FNR > 1 {$NF=barcodes[FNR]; print}' barcoding_summary_sorted.txt sequencing_summary_sorted.txt > dorado_temp
echo "$header" | cat - dorado_temp > dorado_temp2
awk '{$1=$1}1' OFS='\t' dorado_temp2 > dorado_summary.tsv
rm dorado_temp*

# run pycoQC
pycoQC \
--summary_file "${FOLDER}"/Dorado_Basecalling/"${RUN}"/demultiplex/dorado_summary.tsv \
--report_title "${PROJECT}"_"${RUN}"_"${MODEL}"_minqscore"${QSCOREDORADO}"_report \
--html_outfile "${FOLDER}"/QC_Stats/"${RUN}"/PycoQC/"${PROJECT}"_"${RUN}"_run"${RUNDATE}".html
# ---------------------------------

# ---------------------------------
# Cleaning Working Directory
# ---------------------------------
rm -R "${LOCAL_WORK_DIR}"/*
# ---------------------------------

# ---------------------------------
mv "${FOLDER}"/Scripts/P01_DoradoGPU.* "${FOLDER}"/Scripts/Sbatch_log/
# ---------------------------------
