#!/bin/bash
#SBATCH -p workq
#SBATCH -t 5:00:00

#Load binaries
module load bioinfo/Nextflow-v19.04.0

nextflow run polishingPipeline_Main.nf --longReads data/SRR5065181_subreads.fastq.gz --lrPolish racon --lrNum 1 --shortReads data/SRR1706186.fastq.gz --srPolish pilon --srNum 1 --assembly data/assemblage.raw.fa --lineage data/enterobacteriales_odb9 --species E_coli_K12  --outdir ./TestPipeline/ 
