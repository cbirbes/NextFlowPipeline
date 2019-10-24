#!/bin/bash
#SBATCH -p workq
#SBATCH -t 5:00:00

#Load binaries
module load bioinfo/Nextflow-v19.04.0

nextflow run polishingPipeline_Main.nf --longReads data/LongReads/SRR5065181.fastq.gz --lrPolish racon --lrNum 1 --shortReads data/ShortReads/SRR1706186.fastq.gz --srPolish pilon --srNum 1 --srAligner bwa --assembly data/Assembly/assembly.raw.fa --lineage data/enterobacteriales_odb9 --species E_coli_K12  --allSteps true --kat true --outdir ./TestPipeline/ 
