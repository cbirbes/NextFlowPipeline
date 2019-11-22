## Genome assembly polishing pipeline
This nextflow pipeline is a genome assembly polishing pipeline available for short and/or long reads.
It use several different polishing and mapping tools and allow the complete polishing of a genome with the software of your choice.
To use this pipeline, you need at least an initial assembly and short and/or long reads.

## Getting Started

```
git clone --recursive https://github.com/Clement-BIRBES/NextFlowPipeline.git
cd NextFlowPipeline
cd test
```

## Tips to run
The pipeline is optimized to run on the genologin cluster (INRA Toulouse).
If you run it locally or on other cluster, you may need to change some part of the code
To test the pipeline, first make sure that <a
href="https://www.nextflow.io/docs/latest/getstarted.html"> Nextflow is installed </a>
then :
```
nextflow run ../main.nf --longReads data/LongReads/SRR5065181.fastq.gz  --lrPolish racon --lrNum 1 --shortReadsBwa1 data/ShortReads/SRR1706185.fastq.gz --srPolish pilon --srNum 1 --assembly data/Assembly/assembly.raw.fa --outdir ./TestPipeline/ --lineage enterobacteriales_odb9 --species E_coli_K12 --kat true
```

If you're running pipeline on clusters, verify the cluster nextflow version and the version in runTest.sh and :
```
sbatch runTest.sh
```

If you have multiple longreads files you must first concatenate them and pass the final file as argument.



The default value of the pipeline are the one giving the best results on our training data (E.coli K12 with LongReads wtdbg2 assembly).

Long reads used: https://www.ebi.ac.uk/ena/data/view/SRR5065181

Short reads used: https://www.ebi.ac.uk/ena/data/view/SRR1706186

## Parameters
```
Usage:

  The typical command for running the pipeline is as follows:
  nextflow run main.nf --longReads longReads.fq.gz --shortReadsLongranger /path/to/Data/ --assembly assembly.fa [options]


Mandatory arguments (minimum 1 reads file, can't use --shortReadsBwa and --shortReadsLongranger in the same run):
  --longReads       Path to long reads fasta or fastq .gz file. If you have more than one file, concatenate them

  --shortReadsBwa1  Path to short reads file 1, using bwa as aligner
  --shortReadsBwa2  Path to short reads file 2, using bwa as aligner
  --shortReadsLongranger    Path to short reads directory, using longranger as aligner

  --assembly        Path to fasta file of genome assembly to polish


Options:
Long Reads options :
  --lrPolish  			Polisher to use with long reads: wtdbg2, racon, pilon, medaka (default value: racon)

  --lrNum	    			Number of long reads polishing to run (default value: 2)

  --alignerExt      Alignment file extension for polishing with racon: paf (reduce polished output quality, improve speed) or sam (default value: sam)

Short Reads options :
  --srPolish	   		Polisher to use with short reads: pilon, racon, freebayes or wtdbg2 (default value: pilon)

  --srNum	    			Number of short reads polishing to run (default value: 2)

  --chunck          Size of targets chuncks for parallelization with pilon (default value : 10000000)

  --noChanges		   	If "true", polish with pilon until there are no more changes (defaults value: false). Overwrite all short reads options

  --sample          If using "--shortReadsLongranger and a directory containing more than 2 reads files, precise the prefix of reads to analyse (--sample option from longranger align)

  --poolseqSize     Size of the poolSeq sample for poolSeq analysis

  --pattern         Maximum 0 for freebayes poolseq pattern (Eg: 3 = 0/0/0/1/..../1), for poolSeq analysis

Analysis and metrics options :
  --lineage		     	Lineage dataset used for BUSCO (if set : Run Busco quality after the last short reads polishing step)

  --species		     	Reference species to built Augustus annotation during BUSCO steps (default value: generic)

  --allSteps        Run BUSCO after each polishing steps (default value : false)

  --kat             If "true" run kat comp at the end of the pipeline (default value: false) for both long reads and short reads.

Output option :
  --outdir		    	The output directory where the results will be saved (default value : ./results/)
```
