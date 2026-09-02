# Hybrid-Genome-Plasmids-Bacteria
**Get whole chromosome and plasmids Bacteria assemblies from Illumina and Nanopore sequencing data**

## Pipeline steps:  
- ONT Basecalling (***Dorado***, ***PycoQC***)
- Reads Cleaning and QC (***fastp***, ***porechop_ABI***, ***NanoQ***)
- Check for contamination using taxonomic profile (***Kraken2***)
- *DeNovo* assembly of bacterial plasmids (***Plassembler***)
- *DeNovo* assembly of bacterial genome (***Hybracter***)
- Control Quality of assemblies (***QUAST***, ***BUSCO***, ***Bandage***)
- Chromosome and plasmids annotation (***Bakta***)
- Detection and typing of plasmids (***staramr***)
- Antimicrobial Resistance Gene Finder in plasmids (***AMRFinderPlus***)

## Instruction:
- Run the pipeline on Austal HPC (***Slurm***)
- Indicate variable's value in the configuration file ``Config_Pipeline_HybridAssembly_Bacteria_Template.sh`` before submitting job
- Example input file: ``Sample_HybridPipeline_Template.xlsx``

## Main tools used:

1. [https://github.com/gbouras13/plassembler] (https://plassembler.readthedocs.io/en/latest/)

     George Bouras, Anna E. Sheppard, Vijini Mallawaarachchi, Sarah Vreugde, Plassembler: an automated bacterial plasmid assembly tool, Bioinformatics, Volume 39, Issue 7, July 2023, btad409, https://doi.org/10.1093/bioinformatics/btad409.

2. [https://github.com/gbouras13/hybracter] (https://hybracter.readthedocs.io/en/latest/)

    George Bouras, Ghais Houtak, Ryan R Wick, Vijini Mallawaarachchi, Michael J. Roach, Bhavya Papudeshi, Louise M Judd, Anna E Sheppard, Robert A Edwards, Sarah Vreugde - Hybracter: Enabling Scalable, Automated, Complete and Accurate Bacterial Genome    Assemblies. (2024) Microbial Genomics doi: https://doi.org/10.1099/mgen.0.001244.

For information on bacteria genome assembly with Nanopore Data (https://rrwick.github.io/)

# Pipeline figure #
soon...
