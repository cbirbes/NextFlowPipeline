## Genome assembly polishing pipeline
This nextflow pipeline is a genome assembly polishing pipeline available for short and/or long reads.
It use several different polishing and mapping tools and allow the complete polishing of a genome with the software of your choice.
To use this pipeline, you need at least an initial assembly and short and/or long reads.

## Getting Started

```
git clone --recursive https://github.com/Clement-BIRBES/NextFlowPipeline.git
cd NextFlowPipeline
```

## Tips to run
To test the pipeline, first make sure that <a
href="https://www.nextflow.io/docs/latest/getstarted.html"> Nextflow is installed </a>
then :
```
nextflow run main.nf --longReads data/LongReads/SRR5065181.fastq.gz --lrPolish racon --lrNum 1 --shortReads data/ShortReads/SRR1706186.fastq.gz --srPolish pilon --srNum 1 --srAligner bwa --assembly data/Assembly/assembly.raw.fa --lineage data/enterobacteriales_odb9 --species E_coli_K12  --allSteps true --kat true --outdir ./TestPipeline/
```

If you're running pipeline on clusters, verify the cluster nextflow version and the version in runTest.sh and :
```
sbatch runTest.sh
```

To use the pipeline with your own data, it's better to use path to Long reads FILE and path to Short reads DIRECTORY.

The default value of the pipeline are the one giving the best results on our training data (E.coli K12 with LongReads wtdbg2 assembly).

Long reads used: https://www.ebi.ac.uk/ena/data/view/SRR5065181

Short reads used: https://www.ebi.ac.uk/ena/data/view/SRR1706186

## Parameters
```
Usage:
  The typical command for running the pipeline is as follows:
  nextflow run polishingPipeline_Main.nf --longReads 'longReads.fq.gz' --shortReads '/path/to/DemultiplexData/' --assembly 'assembly.fa' [options]

Mandatory arguments:
  --longReads       Path to long reads fasta or fastq .gz file
       AND/OR
  --shortReads      Path to short reads directory (For best performance, prefer 2 reads files/directory)

  --assembly        Fasta file of genome assembly to polish

Options:
Long Reads options :
  --lrPolish        Polisher to use with long reads: wtdbg2, racon, pilon, medaka (default value: racon)

  --lrNum           Number of long reads polishing to run (default value: 2)

  --alignerExt      Alignment file extension for polishing with racon: paf (reduce quality, improve speed) or sam (default value: sam)

Short Reads options :
  --srPolish        Polisher to use with short reads: pilon, racon, freebayes or wtdbg2 (default value: pilon)

  --srNum           Number of short reads polishing to run (default value: 2)

  --chunck          Size of targets for parallelization with pilon (default value : 10000000)

  --noChanges       If "true", polish until there are no more changes (defaults value: false). Overwrite all short reads options

  --srAligner       Aligner to use for short reads: bwa or longranger (default value : longranger). For bwa precise the entire path
                    to the reads, for samtools precise the directory containing the reads.

  --sample          If using "--srAligner longranger" and a directory containing more than 2 reads files, precise the prefix of reads to analyse
                    (--sample option from longranger align)

  --poolseqSize     Size of the poolSeq sample for poolSeq analysis

  --pattern         Maximum 0 for freebayes poolseq pattern (Eg: 3 = 0/0/0/1/..../1), for poolSeq analysis

Analysis and metrics options :
  --lineage         Lineage dataset used for BUSCO (if set : Run Busco quality after the last short reads polishing step)

  --species         Reference species to built Augustus annotation during BUSCO steps (default value: generic)

  --allSteps        Run BUSCO after each polishing steps (default value : false)

  --kat             If "true" run kat comp at the end of the pipeline (default value: false) for both long reads and short reads.

Output option :
  --outdir          The output directory where the results will be saved (default value : ./results/)
```
