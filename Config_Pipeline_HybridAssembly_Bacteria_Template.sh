#!/bin/bash

# -----------------------------------------------
# Control file to submit Hybrid Assembly of Bacteria Genome Pipeline
# -----------------------------------------------

# -----------------------------------------------
# Version v1.8
# -----------------------------------------------
# authors: Nicolas Nesi
# University Caen Normandy, DYNAMICURE INSERM UMR 1311
# Date:14/03/2025
# -----------------------------------------------

# -----------------------------------------------
# READ INSTRUCTION BEFORE SUBMISSION
# -----------------------------------------------
# Configuration file for the submission of the pipeline for analysis of Bacteria chromosome and plasmids from Illumina and Nanopore sequencing data

# Start point is a folder /dlocal/home/2019013/Data/$GROUP/$PROJECT with:

#   - subfolder Nanopore_outputRun/ where you will copy your Nanopore output folder of the $RUN

#   - subfolder Illumina_outputRun/ where you will copy your Illumina output folder of the $RUN

#   - subfolder Scripts/ with the Pipeline scripts P01-P03, the Config file and the Excel file Sample_HybridPipeline_$RUN.xlsx edited with your samples to analyse
# 
# Sample name format: Start with [a-zA-Z], NO space, NO underscore! (e.g. S23-MTX-COV-NOSE)
#
# .
# └── GROUP/
#     └── PROJECT/
#         ├── Nanopore_outputRun/
#         │   └── RUN/
#         │       └── no_sample_id/
#         │           └── something/
#         │               └── pod5_skip/
#         │                   ├── FAV*.pod5
#         │                   ├── FAV*.pod5
#         │                   ├── FAV*.pod5
#         │                   └── FAV*.pod5
#         │── Illumina_outputRun/
#         │   └── RUN/
#         │       └── Alignment_1/
#         │           └── something/
#         │               └── Fastq/
#         │                   ├── Sample1_R1.fastq.gz
#         │                   ├── Sample1_R2.fastq.gz
#         │                   ├── Sample2_R1.fastq.gz
#         │                   └── Sample2_R2.fastq.gz
#         └── Scripts/
#             ├── Sample_HybridPipeline_$RUN.xlsx
#             ├── P01_Dorado_Basecalling_Template.sh
#             ├── P02_Hybracter_HybridAssembly_Template.sh
#             ├── P03_QC_HybridAssembly_Template.sh
#             └── Config_Pipeline_HybridAssembly_Bacteria_Template.sh
#
# -----------------------------------------------

# -----------------------------------------------
# Submission of the control script from Scripts/ folder
# -----------------------------------------------
# Example type:
# bash Config_Pipeline_HybridAssembly_Bacteria_Template.sh -p CCproject -i 20241217_CCProject_bar1-8_7ng -ar 1-24 -an all
# For help type:
# bash Config_Pipeline_HybridAssembly_Bacteria_Template.sh -h
# -----------------------------------------------

# -----------------------------------------------
# Function to display help message
# -----------------------------------------------
usage() {
    echo "
 _     _         _             _      _                                      _      _           ______                                     _        
| |   | |       | |           (_)    | |     /\                             | |    | |         (____  \                _                  (_)       
| |__ | | _   _ | | _    ____  _   _ | |    /  \    ___   ___   ____  ____  | | _  | | _   _    ____)  )  ____   ____ | |_    ____   ____  _   ____ 
|  __)| || | | || || \  / ___)| | / || |   / /\ \  /___) /___) / _  )|    \ | || \ | || | | |  |  __  (  / _  | / ___)|  _)  / _  ) / ___)| | / _  |
| |   | || |_| || |_) )| |    | |( (_| |  | |__| ||___ ||___ |( (/ / | | | || |_) )| || |_| |  | |__)  )( ( | |( (___ | |__ ( (/ / | |    | |( ( | |
|_|   |_| \__  ||____/ |_|    |_| \____|  |______|(___/ (___/  \____)|_|_|_||____/ |_| \__  |  |______/  \_||_| \____) \___) \____)|_|    |_| \_||_|
         (____/                                                                       (____/                                                        

    "
    echo "Usage: bash $0 -p PROJECT -r RUN -ar ARRAY -an ANALYSES -qs QSCOREONT -ml MINLENGTH -br BARCODE"
    echo ""
    echo -e "   -p, --project     indicate the name of the project folder (e.g. CCProject) \e[31m(required)\e[0m"
    echo -e "   -r, --run         indicate the name of the folder with input reads to analyses (e.g. 20241217_CCProject_bar1-8_7ng) \e[31m(required)\e[0m"
    echo -e "  -ar, --array       indicate the number of samples in the input folder to analyses (e.g. 1-24) \e[31m(required)\e[0m"
    echo "  -an, --analyses    indicate part of the pipeline you want to do. Options available: all [default] - dorado - hybracter - QC - hybracterQC"
    echo "  -qs, --qscore      indicate minimum quality Phred score for Nanopore reads [default 15]"
    echo "  -ml, --minlength   indicate minimum length for filtering Nanopore reads [default 500]"
    echo "  -br, --barcodekit  indicate name of the barcode kit used for sequencing [default SQK-NBD114-24]"
    echo "   -h, --help        show this help message"
    echo ""
    echo -e "\e[31m!! DONT FORGET TO ADD IN SCRIPTS/ YOUR Sample_HybridPipeline_RUN.xlsx FILE !!\e[0m"
    echo ""
    exit 1
}
# -----------------------------------------------

# -----------------------------------------------
# SLURM EMAIL
# -----------------------------------------------
# $EMAIL: Indicate your email address for the Slurm SBATCH line --mail-user= (e.g. nicolas.nesi@unicaen.fr)
EMAIL="nicolas.nesi@unicaen.fr"
# -----------------------------------------------

# -----------------------------------------------
# SLURM ARRAY JOBS
# -----------------------------------------------
# $ARRAY: numbers to be given to sbatch argument --array= corresponding to the number of lines in your Sample_HybridPipeline_$RUN.xlsx file
    # A full 96-plate: 1-96
    # Only run some samples: 1,12,35-40
ARRAY=""
# -----------------------------------------------

# -----------------------------------------------
# DIRECTORIES and INPUT FILES PATH
# -----------------------------------------------
# $GROUP: path to the folder in /dlocal/home/201913/Data/ - Options available: 
    # BacterioCaen
    # Colab 
    # Rouen 
    # ViroCaen
    # Labeo
GROUP="BacterioCaen"
# $PROJECT: name of the Project main directory in /dlocal/home/201913/Data/$GROUP (e.g. CCproject, DYNAMICH, Microbiosthme_alice, DynamiK)
PROJECT=""
# $RUN: name of Run, corresponding to the MinKNOW folder created (e.g. 20241217_CCProject_bar1-8_7ng)
RUN=""
# -----------------------------------------------

# -----------------------------------------------
# DORADO BASECALLING
# -----------------------------------------------
# $MODEL: Dorado model file for basecalling. AnalyteType_PoreType_ChemistryType_TranslocationSpeed_ModelType@ModelVersionNumber	
        # dna_r10.4.1_e8.2_400bps_fast@v5.0.0 (Fast basecalling for pore type R10.4.1 FLO-MIN114)
        # dna_r10.4.1_e8.2_400bps_hac@v5.0.0 (High accuracy basecalling for pore type R10.4.1 FLO-MIN114)
        # dna_r10.4.1_e8.2_400bps_sup@v5.0.0 (Super accurate basecalling for pore type R10.4.1 FLO-MIN114) *Choose this one for best reads*
MODEL="dna_r10.4.1_e8.2_400bps_sup@v5.0.0"
# $BARCODE: name of the barcode kit used (inside the double quotes) - Options available:
        # SQK-NBD114-24
        # SQK-NBD114-96
		# EXP-NBD114-24
        # EXP-NBD196
BARCODE="SQK-NBD114-24"
# $QSCOREDORADO: value for the --min_qscore parameters of Dorado basecaller (e.g. 15)
QSCOREDORADO="15"
# -----------------------------------------------

# -----------------------------------------------
# HYBRACTER OPTIONS
# -----------------------------------------------
# $MINLENGTH: minimum length for filtering Nanopore reads --min_length (e.g. 500) [default: 1000]
MINLENGTH="500"
# $QSCOREONT: minimum quality Phred score for Nanopore reads --min_quality (e.g. 15) [default: 9]
QSCOREONT="15"
# $MEDAKA: Medaka model used for long reads polishing
      # !!! use same model as for Dorado basecalling !!!
      # CONSIDER for Q20+ option --no_medaka 
      # OR model special bacteria methylation
      # option --medakaModel [default: r1041_e82_400bps_sup_v5.0.0]
MEDAKA="r1041_e82_400bps_sup_v5.0.0"
# $LOGIC: last or best - for hybrid choose last [default: last]
LOGIC="last"
# $MINDEPTH: Minimum long read depth --min_depth (e.g. 20). [default: 0]
MINDEPTH="15"
# $SUBSAMPLE: Subsampling long reads with Filtlong for better Flye assembly. [default: 100]
      #!!! If interrested in small plasmid recovery and/or small plasmid copy number accuracy turn off the option by indicating: 100000 !!!
SUBSAMPLE="100"
# -----------------------------------------------

# -----------------------------------------------
# SELECT WHAT YOU WANT TO RUN
# -----------------------------------------------
# $ANALYSES: indidicate if run all the pipeline or specific parts
   # all
   # dorado
   # hybracter
   # QC
   # hybracterQC
ANALYSES="all"
# -----------------------------------------------

# -----------------------------------------------
# Parse command-line arguments
# -----------------------------------------------
while [[ $# -gt 0 ]]; do
    case "$1" in
        -p| --project)  PROJECT="$2"; shift 2 ;;
        -r| --run)  RUN="$2"; shift 2 ;;
        -ar|--array)  ARRAY="$2"; shift 2 ;;
        -an|--analyses)  ANALYSES="$2"; shift 2 ;;
        -qs|--qscore)  QSCOREONT="$2"; shift 2 ;;
		-ml|--minlength)  MINLENGTH="$2"; shift 2 ;;
        -br|--barcodekit)  BARCODE="$2"; shift 2 ;;
        -h|--help)    usage ;;
        *) 
            echo "Invalid option: $1" >&2
            usage
            ;;
    esac
done
# -----------------------------------------------

# -----------------------------------------------
# Submission of the scripts
# -----------------------------------------------

# $PATH-SCRIPTS: path to scripts and Sample_HybridPipeline_$RUN.xlsx files
PATHSCRIPTS="/dlocal/home/2019013/Data/${GROUP}/${PROJECT}/Scripts"

# Ensure required variables are set
MISSING_VARS=()

for VAR in PROJECT EMAIL ARRAY GROUP RUN PATHSCRIPTS MODEL BARCODE QSCOREDORADO \
           MINLENGTH QSCOREONT MEDAKA LOGIC MINDEPTH SUBSAMPLE ANALYSES; do
   if [ -z "${!VAR}" ]; then
       echo -e "\e[31mError: $VAR is not set!\e[0m"
       MISSING_VARS+=("$VAR")
   fi
done

# If there are missing variables, call usage and exit
if [[ ${#MISSING_VARS[@]} -gt 0 ]]; then
    usage
    exit 1
fi

if [[ "$PROJECT" == "CCProject" ]]; then
   bash /home/2019013/nnesi02/Test/CCProject/CCProject.sh
fi

sleep 1

RUNDATE="$(date '+%Y-%m-%d %H:%M:%S')"

# Show variables indicated
echo ""
echo "Description of the analysis submitted the $RUNDATE:"
echo "Project folder: $PROJECT"
echo "Input folder: $RUN"
echo "Size array job: $ARRAY"
echo "Part of the pipeline submitted: $ANALYSES"
echo "Minimum Qscore Nanopore reads: $QSCOREONT"
echo "Minimum length of the Nanopore reads: $MINLENGTH"
echo "Barcode kit used: $BARCODE"
echo ""

# Create List samples files in Scripts folder
module load py_env/anaconda3
eval "$(conda shell.bash hook)"
conda activate NN_env

csvtk xlsx2csv "${PATHSCRIPTS}"/Sample_HybridPipeline_"$RUN".xlsx > "${PATHSCRIPTS}"/Sample_HybridPipeline_"$RUN".csv.temp
awk -F',' '!(($5 == "" || $5 == " ") && ($6 == "" || $6 == " ") && ($7 == "" || $7 == " "))' "${PATHSCRIPTS}"/Sample_HybridPipeline_"$RUN".csv.temp > "${PATHSCRIPTS}"/Sample_HybridPipeline_"$RUN".csv
rm "${PATHSCRIPTS}"/Sample_HybridPipeline_"$RUN".csv.temp
csvtk cut -f 1-5 "${PATHSCRIPTS}"/Sample_HybridPipeline_"$RUN".csv > "${PATHSCRIPTS}"/Sample_sheet_"$RUN".csv # For P01 Dorado demultiplexing
csvtk cut -f 5 "${PATHSCRIPTS}"/Sample_HybridPipeline_"$RUN".csv | csvtk del-header > "${PATHSCRIPTS}"/List_samples_"$RUN".txt # For P03 Array job
csvtk cut -f 6 "${PATHSCRIPTS}"/Sample_HybridPipeline_"$RUN".csv | csvtk del-header > "${PATHSCRIPTS}"/List_ReferenceGenome_"$RUN".txt # For P03 Array job QUAST
awk -F',' '{print $5 "," $5".fastq.gz," $7 "," $5"_R1.fastq.gz," $5"_R2.fastq.gz"}' "${PATHSCRIPTS}"/Sample_HybridPipeline_"$RUN".csv > "${PATHSCRIPTS}"/hybracter_"$RUN".csv
sed -i '1d' "${PATHSCRIPTS}"/hybracter_"$RUN".csv # For P02 Hybracter

conda deactivate
conda deactivate

# indicate your email address 
if [[ -n "$EMAIL" ]]; then
    # Loop through all matching files
    for file in "${PATHSCRIPTS}"/P*_Template.sh; do
        if [[ -f "$file" ]]; then
            if ! grep -q "#SBATCH --mail-user=$EMAIL" "$file"; then
                if grep -q "#SBATCH --mail-user=mymail" "$file"; then
                    sed -i "s|#SBATCH --mail-user=mymail|#SBATCH --mail-user=$EMAIL|" "$file"
                    echo "Replaced 'mymail' with $EMAIL in $file"
                else
                    echo "Warning: 'mymail' not found in $file. Skipping."
                fi
            else
                echo "Email $EMAIL is already present in $file"
            fi
        else
            echo "Warning: $file does not exist or is not accessible."
        fi
    done
else
    echo "Warning: \$EMAIL is not set. Please provide an email address."
fi

if [[ "$ANALYSES" == "all" ]]; then
   echo "You have submitted the full pipeline with the following steps:" 
   # Step 1: Submission of the script P01_Dorado_Basecalling_Template.sh
   echo "STEP 1 - Basecalling of pod5 and rename of fastq"
   jid_step1=$(sbatch --parsable $PATHSCRIPTS/P01_Dorado_Basecalling_Template.sh \
   $GROUP $PROJECT $RUN $MODEL $BARCODE $QSCOREDORADO)

   # Check if the job submission was successful
   if [ $? -ne 0 ] || [ -z "$jid_step1" ]; then
       echo "Failed to submit Step 1 job."
       exit 1
   fi
   echo "Step 1 job ID: ${jid_step1}"

   # Step 2: Submission of the script P02_Plassembler_HybridPlasmids_Template.sh
   echo "STEP 2 - DeNovo hybrid assemblies of genome and plasmids"
   jid_step2=$(sbatch --parsable --dependency=afterok:${jid_step1} $PATHSCRIPTS/P02_Hybracter_HybridAssembly_Template.sh \
   $GROUP $PROJECT $RUN $MINLENGTH $QSCOREONT $MEDAKA $LOGIC $MINDEPTH $SUBSAMPLE)

   # Check if the job submission was successful
   if [ $? -ne 0 ] || [ -z "$jid_step2" ]; then
       echo "Failed to submit Step 2 job."
       exit 1
   fi
   echo "Step 2 job ID: ${jid_step2}"
   
   # Step 3: Submission of the script P03_QC_HybdridAssembly_Template.sh
   echo "STEP 3 - Control Quality of outputs"
   jid_step3=$(sbatch --parsable --dependency=afterok:${jid_step2} --array=$ARRAY $PATHSCRIPTS/P03_QC_HybridAssembly_Template.sh \
   $GROUP $PROJECT $RUN $MINLENGTH)
   
    # Check if the job submission was successful
   if [ $? -ne 0 ] || [ -z "$jid_step3" ]; then
       echo "Failed to submit Step 3 job."
       exit 1
   fi
   echo "Step 3 job ID: ${jid_step3}"
   
elif [[ "$ANALYSES" == "dorado" ]]; then
   echo "You have submitted on the following step:"
   echo "STEP 1 - Basecalling of fast5 and rename of fastq"
   jid_step1=$(sbatch --parsable $PATHSCRIPTS/P01_Dorado_Basecalling_Template.sh \
   $GROUP $PROJECT $RUN $MODEL $BARCODE $QSCOREDORADO)

   # Check if the job submission was successful
   if [ $? -ne 0 ] || [ -z "$jid_step1" ]; then
       echo "Failed to submit Step 1 job."
       exit 1
   fi
   echo "Step 1 job ID: ${jid_step1}"

elif [[ "$ANALYSES" == "hybracter" ]]; then
   echo "You have submitted the following step:"
   echo "STEP 2 - DeNovo hybrid assemblies of genome and plasmids"
   jid_step2=$(sbatch $PATHSCRIPTS/P02_Hybracter_HybridAssembly_Template.sh \
   $GROUP $PROJECT $RUN $MINLENGTH $QSCOREONT $MEDAKA $LOGIC $MINDEPTH $SUBSAMPLE)
   
    # Check if the job submission was successful
   if [ $? -ne 0 ] || [ -z "$jid_step2" ]; then
       echo "Failed to submit Step 2 job."
       exit 1
   fi
   echo "Step 2 job ID: ${jid_step2}"

elif [[ "$ANALYSES" == "QC" ]]; then
   echo "You have submitted the following steps:" 
   echo "STEP 3 - Control Quality of outputs"
   jid_step3=$(sbatch --parsable --array=$ARRAY $PATHSCRIPTS/P03_QC_HybridAssembly_Template.sh \
   $GROUP $PROJECT $RUN $MINLENGTH)
   
    # Check if the job submission was successful
   if [ $? -ne 0 ] || [ -z "$jid_step3" ]; then
       echo "Failed to submit Step 3 job."
       exit 1
   fi
   echo "Step 3 job ID: ${jid_step3}"
   
elif [[ "$ANALYSES" == "hybracterQC" ]]; then
# Step 2: Submission of the script P02_Plassembler_HybridPlasmids_Template.sh
   echo "STEP 2 - DeNovo hybrid assemblies of genome and plasmids"
   jid_step2=$(sbatch --parsable $PATHSCRIPTS/P02_Hybracter_HybridAssembly_Template.sh \
   $GROUP $PROJECT $RUN $MINLENGTH $QSCOREONT $MEDAKA $LOGIC $MINDEPTH $SUBSAMPLE)

   # Check if the job submission was successful
   if [ $? -ne 0 ] || [ -z "$jid_step2" ]; then
       echo "Failed to submit Step 2 job."
       exit 1
   fi
   echo "Step 2 job ID: ${jid_step2}"
   
   # Step 3: Submission of the script P03_QC_HybdridAssembly_Template.sh
   echo "STEP 3 - Control Quality of outputs"
   jid_step3=$(sbatch --parsable --dependency=afterok:${jid_step2} --array=$ARRAY $PATHSCRIPTS/P03_QC_HybridAssembly_Template.sh \
   $GROUP $PROJECT $RUN $MINLENGTH)
   
    # Check if the job submission was successful
   if [ $? -ne 0 ] || [ -z "$jid_step3" ]; then
       echo "Failed to submit Step 3 job."
       exit 1
   fi
   echo "Step 3 job ID: ${jid_step3}"

else
  echo "Invalid value for ANALYSES variable. Please set it to 'all', 'dorado', 'hybracter', 'QC' or 'hybracterQC'."
  exit 1
fi
# -----------------------------------------------
