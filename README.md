# Polishing Pipeline
This nextflow pipeline is an assembly polishing pipeline available for short and / or long reads.
It use several tools like pilon, racon, medaka, freebayes or wtdbg2.
To use it, you need at least an initial assembly and short and / or long reads.

## Parameters
Usage:
The typical command for running the pipeline is as follows:
nextflow run polishingPipeline_Main.nf --longReads 'longReads.fq.gz' --shortReads '/path/to/DemultiplexData/' --assembly 'assembly.fa' [options]

Mandatory arguments:
--longReads       Path to long reads fasta or fastq .gz file
     AND/OR
--shortReads      Path to demultiplexed 10X short reads directory

--assembly        Fasta file of genome assembly to polish

Options:

--LRPolish  			Polisher to use with long reads: wtdbg2, racon, pilon, medaka (default value: racon)

--LRNum	    			Number of long reads polishing to run (default value: 2)

--AlignerExt      Alignment file extension for polishing with racon: paf (reduce quality, improve speed) or sam (default value: sam)

--SRPolish	   		Polisher to use with short reads: pilon, racon, freebayes or wtdbg2 (default value: pilon)

--SRNum	    			Number of short reads polishing to run (default value: 2)

--Chunck          Size of targets for parallelization with pilon (default value : 10000000)

--NoChanges		   	If "true", polish until there are no more changes (defaults value: false). Overwrite all short reads options

--Outdir		    	The output directory where the results will be saved (default value : ./results/)

--Lineage		     	Lineage dataset used for BUSCO (Run Busco quality at the end of the pipeline if set)

--Species		     	Reference species to built Augustus annotation during BUSCO (default value: generic)

--Kat             If "true" run kat comp at the end of the pipeline (default value: false)

--SRAligner			  Aligner to use for short reads: samtools or longranger (default value : longranger)
