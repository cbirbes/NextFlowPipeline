# Polishing Pipeline
This nextflow pipeline is an assembly polishing pipeline available for short and / or long reads.
It use several tools like pilon, racon, medaka, freebayes or wtdbg2.
To use it, you need at least an initial assembly and short and / or long reads.

## Getting Started

git clone https://github.com/Clement-BIRBES/NextFlowPipeline.git
cp NextFlowPipeline/* /path/to/your/working/directory

## Parameters
Usage:
The typical command for running the pipeline is as follows:
nextflow run polishing --longReads 'longReads.fq.gz' --shortReads '/path/to/DemultiplexData/' --assembly 'assembly.fa'

Mandatory arguments:
--longReads       Path to long reads fasta or fastq .gz file
     AND/OR
--shortReads      Path to demultiplexed 10X short reads directory

--assembly        Fasta file of genome assembly to polish

Options:

--LRPolish  			Polisher to use for long reads: wtdbg2, racon, pilon, medaka (default value: 'racon')

--LRNum	    			Specify the number of long reads polishing to run (default value: 2)

--AlignerExt      Chose aligner extension for racon polishing: 'paf' (.paf reduce quality, improve speed) or 'sam' (default value: 'sam')

--SRPolish	   		Polisher to use for short reads: pilon, racon, freebayes or wtdbg2 (default value: 'pilon')

--SRNum	    			Specify the number of short reads polishing to run (default value: 2)

--NoChanges		   	true : Pilon polishing until there are no more assembly changes (defaults value: false). Overwrite SRPolish and SRNum options

--outdir		    	The output directory where the results will be saved (default value : './results/')

--lineage		     	Lineage dataset used for BUSCO (Run Busco quality if set)

--species		     	Reference species to built Augustus annotation during BUSCO (default value: 'generic')

--kat             true : kat evaluation (default value: "false")

--reference			  Reference genome used for Quast comparison

--genes           Gene and operon annotations used for Quast

--chunck          Contig length for pilon Parallelization (default value : 10000000)

--SRAligner			  Aligner to use for short reads: samtools or longranger (default value : longranger)

--poolseqSize     Size of the poolSeq sample

--pattern         Maximum 0 for freebayes poolseq pattern (Eg: 3 = 0/0/0/1/..../1)
