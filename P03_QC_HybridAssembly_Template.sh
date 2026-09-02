#!/bin/bash
#SBATCH -J P03_QC                           # Job name 
#SBATCH -o P03_QC.%a.%A.out                 # File to which stdout will be written
#SBATCH -e P03_QC.%a.%A.err                 # File to which stderr will be written
#SBATCH -p 2tcourt                          # Partition
#SBATCH -N 1                                # Number of node
#SBATCH -c 24                               # Number of cores/cpus
#SBATCH -t 00-01:00                         # Runtime in DD-HH:MM
#SBATCH --mem 400G                          # Memory for all cores in Gbytes
#SBATCH --mail-type=ALL                     # BEGIN,END,ALL
#SBATCH --mail-user=mymail                  # Email address

# ---------------------------------
# Version v1.8
# ---------------------------------
# authors: Nicolas Nesi
# University Caen Normandy, DYNAMICURE INSERM UMR 1311
# Date:30/12/2024
# ---------------------------------

# ---------------------------------
# Variables and paths
# ---------------------------------
# Positional parameters
GROUP="${1}"
PROJECT="${2}"
RUN="${3}"
MINLENGTH="${4}"

# Paths
FOLDER="/dlocal/home/2019013/Data/$GROUP/$PROJECT"
PATH1="$FOLDER/RawReads_Nanopore/$RUN"
PATH2="$FOLDER/RawReads_Illumina/$RUN"
PATH3="$FOLDER/CleanedReads_Nanopore/$RUN"
PATH4="$FOLDER/CleanedReads_Illumina/$RUN"
TARGET="${FOLDER}/QC_Stats/${RUN}/Kraken2/"

# Databases
KRAKENDB="/dlocal/home/2019013/Databases/Kraken2_PlusPFP-2024"
QUASTDB="/dlocal/home/2019013/Databases/QUAST/Reference_genomes"
BUSCODB="/dlocal/home/2019013/Databases/BUSCO/lineages"

# Array job
SAMPLE="$(sed -n -e "${SLURM_ARRAY_TASK_ID}p" $FOLDER/Scripts/List_samples_$RUN.txt)"
REFERENCE="$(sed -n -e "${SLURM_ARRAY_TASK_ID}p" $FOLDER/Scripts/List_ReferenceGenome_$RUN.txt)"
RUNDATE="$(date '+%Y-%m-%d')"

# Threads
export THREADS="${SLURM_CPUS_PER_TASK}"
echo "number of threads used ${THREADS}"
# ---------------------------------

# ---------------------------------
echo Analysis run the: "${RUNDATE}"
echo Working directory: "${PWD}"
echo Submitted with Variables: Group "${GROUP}" Project "${PROJECT}" Run "${RUN}" Minimum length Nanopore reads "${MINLENGTH}"
# ---------------------------------

# ---------------------------------
# QC assembly
# ---------------------------------
module purge
module load py_env/anaconda3
eval "$(conda shell.bash hook)"
conda activate denovo_env
# include:
# quast=5.3.0 

echo "Running QUAST on hybracter assembly for sample: $SAMPLE"

quast \
"${FOLDER}/Hybracter/${RUN}/FINAL_OUTPUT/complete/${SAMPLE}_final.fasta" \
--output-dir "${FOLDER}/QC_Stats/${RUN}/QUAST/${SAMPLE}_QUAST" \
--nanopore "${PATH3}/${SAMPLE}_cleaned_min${MINLENGTH}.fastq.gz" \
--pe1 "${PATH4}/${SAMPLE}_cleaned_R1.fastq.gz" \
--pe2 "${PATH4}/${SAMPLE}_cleaned_R2.fastq.gz" \
-r "${QUASTDB}"/"${REFERENCE}"_genomic.fna \
--features gene:"${QUASTDB}"/"${REFERENCE}"_genomic.gff \
--gene-finding \
--plots-format svg \
--threads "${THREADS}" 2>&1 | tee -a "${FOLDER}/Log/${RUN}_run${RUNDATE}/${SAMPLE}_QUAST_run${RUNDATE}.log"


conda activate BUSCO_env
# Include:
# busco=5.8.0

busco \
--in "${FOLDER}/Chromosome/${RUN}/${SAMPLE}_chromosome.fasta" \
--mode genome \
--offline \
--tar \
--force \
--lineage_dataset "${BUSCODB}"/bacteria_odb10/ \
--out_path "${FOLDER}/QC_Stats/${RUN}/BUSCO/" \
--out "${SAMPLE}"_BUSCO/ \
--cpu "${THREADS}"

python3 /soft/2019013/Conda_env/BUSCO_env/bin/generate_plot.py --working_directory "${FOLDER}/QC_Stats/${RUN}/BUSCO/${SAMPLE}"_BUSCO/
mv "${FOLDER}/QC_Stats/${RUN}/BUSCO/${SAMPLE}"_BUSCO/busco_figure.png "${FOLDER}/QC_Stats/${RUN}/BUSCO/${SAMPLE}"_BUSCO/"${SAMPLE}"_busco_figure.png

rm "${FOLDER}"/Scripts/busco*.log
# ---------------------------------

# ---------------------------------
# Statistic Nanopore reads
# ---------------------------------
conda activate clean-ONT_env
# include:
# nanoq=0.10.0 

echo "Starting statistic of Nanopore raw fastq for sample: $SAMPLE"

nanoq \
-i "${PATH1}/${SAMPLE}.fastq.gz" \
--stats \
--header \
-vvv > "$FOLDER/QC_Stats/$RUN/RawReads_Nanopore/${SAMPLE}_RawReads_Nanopore_Report.txt"

gawk -i inplace 'BEGINFILE {split(FILENAME,f,".")} NR==1{print f[1];print;next} {print}' "$FOLDER/QC_Stats/$RUN/RawReads_Nanopore/${SAMPLE}_RawReads_Nanopore_Report.txt"
sed -i '2,5d' "$FOLDER/QC_Stats/$RUN/RawReads_Nanopore/${SAMPLE}_RawReads_Nanopore_Report.txt"
sed -i 's/_RawReads_Report//' "$FOLDER/QC_Stats/$RUN/RawReads_Nanopore/${SAMPLE}_RawReads_Nanopore_Report.txt"
sed -i 's/ \+ /\t/g' "$FOLDER/QC_Stats/$RUN/RawReads_Nanopore/${SAMPLE}_RawReads_Nanopore_Report.txt"
sed -i "s|${FOLDER}/QC_Stats/${RUN}/RawReads_Nanopore/||g" "$FOLDER/QC_Stats/$RUN/RawReads_Nanopore/${SAMPLE}_RawReads_Nanopore_Report.txt"

echo "Starting statistic of Nanopore cleaned fastq for sample: $SAMPLE"

nanoq \
-i "${PATH3}/${SAMPLE}_cleaned_min${MINLENGTH}.fastq.gz" \
--stats \
--header \
-vvv > "$FOLDER/QC_Stats/$RUN/CleanedReads_Nanopore/${SAMPLE}_cleaned_min${MINLENGTH}_Nanopore_Report.txt"

gawk -i inplace 'BEGINFILE {split(FILENAME,f,".")} NR==1{print f[1];print;next} {print}' "$FOLDER/QC_Stats/$RUN/CleanedReads_Nanopore/${SAMPLE}_cleaned_min${MINLENGTH}_Nanopore_Report.txt"
sed -i '2,5d' "$FOLDER/QC_Stats/$RUN/CleanedReads_Nanopore/${SAMPLE}_cleaned_min${MINLENGTH}_Nanopore_Report.txt"
sed -i "s/^$SAMPLE.*/$SAMPLE/" "$FOLDER/QC_Stats/$RUN/CleanedReads_Nanopore/${SAMPLE}_cleaned_min${MINLENGTH}_Nanopore_Report.txt"
sed -i 's/ \+ /\t/g' "$FOLDER/QC_Stats/$RUN/CleanedReads_Nanopore/${SAMPLE}_cleaned_min${MINLENGTH}_Nanopore_Report.txt"
sed -i "s|${FOLDER}/QC_Stats/${RUN}/CleanedReads_Nanopore/||g" "$FOLDER/QC_Stats/$RUN/CleanedReads_Nanopore/${SAMPLE}_cleaned_min${MINLENGTH}_Nanopore_Report.txt"
# ---------------------------------

# ---------------------------------
# Taxonomic Profile
# ---------------------------------
echo "Starting taxonomic profile for sample $SAMPLE"

conda deactivate
conda activate MTX_env
# include:
# Kraken2=2.1.3
# KrakenTools=1.2
# Krona=2.8.1

kraken2 \
--db "${KRAKENDB}" \
--memory-mapping \
--threads "${THREADS}" \
--minimum-hit-groups 3 \
--output "${FOLDER}/QC_Stats/${RUN}/Kraken2/${SAMPLE}_ONT_Kraken2.tsv" \
--report "${FOLDER}/QC_Stats/${RUN}/Kraken2/${SAMPLE}_ONT.k2report" \
--unclassified-out "${FOLDER}/QC_Stats/${RUN}/Kraken2/${SAMPLE}_ONT_NoHit.fastq.gz" \
"${PATH3}/${SAMPLE}_cleaned_min${MINLENGTH}.fastq.gz" 2>&1 | tee -a "${FOLDER}/Log/${RUN}_run${RUNDATE}/${SAMPLE}_Kraken2_ONT.log"

# Kraken2 Krona
kreport2krona.py \
-r "${FOLDER}/QC_Stats/${RUN}/Kraken2/${SAMPLE}_ONT.k2report" \
-o "${FOLDER}/QC_Stats/${RUN}/Kraken2/${SAMPLE}_ONT_Kraken2_Krona.txt"

ktImportText "${FOLDER}/QC_Stats/${RUN}/Kraken2/${SAMPLE}_ONT_Kraken2_Krona.txt" \
-o "${FOLDER}/QC_Stats/${RUN}/Kraken2/${SAMPLE}_ONT_Kraken2_Krona.html"

# remove empty NoHit fastq.gz
find "$TARGET" -type f -name "*.fastq.gz" -size 0 -exec rm -f {} +
# ---------------------------------

# ---------------------------------
# Bandage
# ---------------------------------
echo "Calculate statistics of assembly for sample: $SAMPLE"

conda deactivate
conda activate bandage_env
# include:
# bandage_ng=2022.09

cp "${FOLDER}/Hybracter/${RUN}/processing/assemblies/${SAMPLE}/assembly_graph.gfa" "${FOLDER}/QC_Stats/${RUN}/gfa/${SAMPLE}_flye_assembly.gfa"
cp "${FOLDER}/Hybracter/${RUN}/processing/plassembler/${SAMPLE}/plassembler_plasmids.gfa" "${FOLDER}/QC_Stats/${RUN}/gfa/${SAMPLE}_plassembler_plasmids.gfa"

Bandage info "${FOLDER}/QC_Stats/${RUN}/gfa/${SAMPLE}_flye_assembly.gfa" > "${FOLDER}/QC_Stats/${RUN}/gfa/${SAMPLE}_flye_assembly_Bandage_statistics.tsv"
Bandage info "${FOLDER}/QC_Stats/${RUN}/gfa/${SAMPLE}_plassembler_plasmids.gfa" > "${FOLDER}/QC_Stats/${RUN}/gfa/${SAMPLE}_plassembler_plasmids_Bandage_statistics.tsv"
# ---------------------------------

# ---------------------------------
mv "${FOLDER}"/Scripts/P03_QC."${SLURM_ARRAY_TASK_ID}".* "${FOLDER}"/Scripts/Sbatch_log
# ---------------------------------
