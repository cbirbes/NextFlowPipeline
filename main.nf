#!/usr/bin/env nextflow
/*
========================================================================================
    			A S S E M B L Y   P O L I S H I N G   P I P E L I N E
========================================================================================
 Assembly Polishing Pipeline. Started October 2018.
 #### Authors ####
 Margot ZAHM <margot.zahm@inra.fr>
 Clément BIRBES <clement.birbes@inra.fr>
 https://github.com/Clement-BIRBES/NextFlowPipeline
----------------------------------------------------------------------------------------
*/
/*
*========================================================
* 						HELP MESSAGE
*========================================================
*/

def helpMessage() {
    log.info"""
	=========================================
	    Assembly polishing Pipeline v1.2
	=========================================
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

    """.stripIndent()
}



/*
*========================================
*=    SET UP CONFIGURATION VARIABLES		=
*========================================
************************
//* Show help message  *
************************
*/
params.help = false
if (params.help){
    helpMessage()
    exit 0
}

//*********************************
//* Default values initialization *
//*********************************
params.longReads		    = false
params.shortReadsBwa1		= false
params.shortReadsBwa2		= false
params.shortReadsLongranger	= false
params.assembly			= false

params.lrPolish 		= 'racon'
params.lrNum			  = 2
LR_number				    = params.lrNum
params.alignerExt   = 'sam'

params.srPolish 		= 'pilon'
params.srNum		    = 2
SR_number				    = params.srNum
params.chunck			  = 10000000
chunckSize          = Channel.value(params.chunck)
params.noChanges		= false
params.sample       = false

params.outdir			  = './results/'

params.lineage			= false
params.species			= 'generic'
params.allSteps     = false

params.kat          = false

params.poolseqSize  = false
PSS                 = params.poolseqSize
params.pattern      = false
Patt                = params.pattern


//******************************************
//* Check if reads are given and not empty *
//******************************************
if (!params.longReads && !params.shortReadsBwa1 && !params.shortReadsBwa2 && !params.shortReadsLongranger) {
	exit 1, "You must specify at least one short or long read file using --shortReadsBwa1, --shortReadsBwa2, --shortReadsLongranger or --longReads."
}

if (params.longReads) {
	LongReads_ch=Channel.fromPath(params.longReads)
	                 .ifEmpty {exit 2, "Long Reads file not found: ${params.longReads}"}
	mode1=true
} else { mode1=false }

if (params.shortReadsLongranger || params.shortReadsBwa1 || params.shortReadsBwa2) {
  if (params.shortReadsLongranger){
	  ShortReads_ch=Channel.fromPath(params.shortReadsLongranger, type: 'dir')
				.ifEmpty {exit 2, "Short Reads file not found: ${params.shortReadsLongranger}"}
    ShortReads_ch2=Channel.value("")
  }
  if (params.shortReadsBwa1){
    ShortReads_ch=Channel.fromPath(params.shortReadsBwa1)
       .ifEmpty {exit 2, "Short Reads file not found: ${params.shortReadsBwa1}"}
    if (params.shortReadsBwa2){
      ShortReads_ch2=Channel.fromPath(params.shortReadsBwa2)
          .ifEmpty {exit 2, "Short Reads file not found: ${params.shortReadsBwa2}"}
    }
    else{
      ShortReads_ch2=Channel.value("")
    }
  }
	mode2=true
} else { mode2=false }


//****************************************************
//* Check if assembly file is given and is not empty *
//****************************************************
if (!params.assembly) {
	exit 1, "You must specify an assembly file."
} else {
	Assembly_ch=Channel.fromPath(params.assembly)
	    			.ifEmpty {exit 1, "Assembly file not found: ${params.assembly}"}
}


//********************************
//* Check BUSCO input parameters *
//********************************
if (params.lineage) {
  Lineage_ch=Channel.fromPath("/usr/local/bioinfo/src/BUSCO/datasets/${params.lineage}", type: 'dir')
						.ifEmpty {exit 1, "Assembly file not found: /usr/local/bioinfo/src/BUSCO/datasets/${params.lineage}"}
            .into {LineageSR_ch; LineageLR_ch; LineageSR1_ch; LineageSR2_ch; LineageSR3_ch}
	BUSCOspecies_ch = Channel.value(params.species)
	mode3 = true
} else {
	mode3 = false
}


//******************************************************
//* Check consistency of Long Reads Polishing argument *
//******************************************************
if (params.lrPolish != 'wtdbg2' && params.lrPolish != 'racon' && params.lrPolish != 'pilon' && params.lrPolish != 'medaka'){
	exit 1, "You must specify an available short reads polisher (wtdbg2, racon, medaka or pilon)"
}
if (params.lrNum.getClass() == Integer) {
	if  (params.lrNum > 5) {
		println "WARN : you are running ${params.lrNum} long reads polishing"
	} else if (params.lrNum <=0) {
		exit 1, "Option --lrNum must be higher than 0"
	}
} else if (params.lrNum.getClass() != Integer) {
	exit 3, "Option --lrNum must be an integer : ${params.lrNum}"
}


//********************************************************
//* Params setting and checking without noChanges option *
//********************************************************
if (!params.noChanges){
//*************************************
//* Check consistency of SR Polishing *
//*************************************
	if (params.srPolish != 'wtdbg2' && params.srPolish != 'racon' && params.srPolish != 'pilon' && params.srPolish != 'freebayes'){
		exit 1, "You must specify an available short reads polisher (wtdbg2, racon, freebayes or pilon)"
	}
	if (params.srNum.getClass() == Integer) {
		if  (params.srNum > 3) {
			println "WARN : you are running ${params.srNum} short reads polishing"
      if (params.srPolish == 'pilon'){
        exit 1, "If you want to polish more than 3 times with pilon, edit the code"
      }
		} else if (params.srNum <=0) {
			exit 1, "Option --srNum must be higher than 0"
		}
	} else if (params.srNum.getClass() != Integer) {
		exit 3, "Option --srNum must be an integer : ${params.srNum}"
	}

//* Sample option for longranger
  if (params.shortReadsLongranger && params.sample) {
    param = '--sample='
    Sample=Channel.value(param.concat(params.sample))
  }
  else {
    Sample=Channel.value('')
  }

  if (params.shortReadsLongranger && params.srPolish == 'racon'){
    exit 1, "Incompatible short reads aligner (longranger) and short reads polisher (racon) combinaison, try bwa as aligner or an other polisher"
  }
}


//**************************
//* Chunck size evaluation *
//**************************
if (params.chunck < 1000000){
	exit 1, "Using small chunck is not recommended, try to run pipeline with a chunk size bigger than 1,000,000"
}
if (params.chunck > 100000000){
	exit 1, "Using big chunck is not recommended, try to run pipeline with a chunk size smaller than 100,000,000"
}


//**************************************
//* Freebayes poolSeq pattern creation *
//**************************************
if (params.pattern != false && params.poolseqSize != false) {
  truePattern = []
  leftSize = PSS.toInteger()-Patt.toInteger()
  for (i=0; i<Patt.toInteger(); i++) (
    truePattern << "0/"
  )
  for (i=0; i<leftSize-1; i++) (
    truePattern << "1/"
  )
  truePattern << "1"
  stringPattern = truePattern.join(",")
}


/*
*========================================================
* 				LONG READS POLISHING PART
*========================================================
//**********************
//* CONFIGURATION PART *
//**********************
* If Long reads are given *
*/
if (mode1) {
//* Create iteration condition, PolisherAssembly channel and the loop for multiple polish using long reads *
	iteration_polisherLR=[]
	condition = { iteration_polisherLR.size()<LR_number ? it : Channel.STOP }

	PolisherAssembly_ch = Channel.create()

	Assembly_ch.mix(PolisherAssembly_ch.map(condition))
			.into{AssemblyMinimap_ch; AssemblyPolisher_ch}


//****************
//* PROCESS PART *
//****************
//* Rename long reads file and add it to channels *
	process rename_long_reads {
		publishDir "${params.outdir}"

		input:
		file longreads from LongReads_ch.collect()

		output:
		file 'longreads.rename.fq.gz' into LongReadsMinimap_ch, LongReadsPolisher_ch, LongReadsKat_ch

		shell:
		if ("${longreads}".endsWith('.gz')){
			'''
			zcat !{longreads} | awk 'BEGIN {i=1} NR%4==1 {print "@read_"i; i++} NR%4!=1 {print $0}' | gzip > longreads.rename.fq.gz
			'''
		} else {
			'''
			cat !{longreads} | awk 'BEGIN {i=1} NR%4==1 {print "@read_"i; i++} NR%4!=1 {print $0}' | gzip > longreads.rename.fq.gz
			'''
		}
	}

//* Map long reads to assembly using Minimap2 with .sam output extension *
	if (params.lrPolish != "racon" || params.alignerExt == "sam" ){
		process minimapLR {
			label 'minimap'

			input:
			file reads from LongReadsMinimap_ch
			file assembly from AssemblyMinimap_ch

			output:
			file 'map.polisher.sam' into MapMinimap_ch

			script:
			"""
      module load bioinfo/minimap2-2.11
			minimap2 -a -t 8 -x map-ont ${assembly} ${reads} > map.polisher.sam
			"""
		}
	}

//* Map long reads to assembly using Minimap2 with .paf output extension (Racon only) *
  if (params.lrPolish == "racon" && params.alignerExt == "paf" ) {
    process minimapRac {
      label 'minimap'

      input:
      file reads from LongReadsMinimap_ch
      file assembly from AssemblyMinimap_ch

      output:
      file 'map.polisher.paf' into MapMinimap_ch

      script:
      """
      module load bioinfo/minimap2-2.11
      minimap2 -t 8 -x map-ont ${assembly} ${reads} > map.polisher.paf
      """
    }
  }

//* Convert .sam file to .sorted.bam file for pilon and wtdbg2
  if (params.lrPolish == "wtdbg2" || params.lrPolish == "pilon"){
    process samtoolsLR {
      label 'samtools'
      input:
      file map from MapMinimap_ch

      output:
      set file("map.polisher.sorted.bam"), file("map.polisher.sorted.bam.bai") into MapSamtoolsLR_ch

      script:
      """
      module load bioinfo/samtools-1.4
      samtools view -S -b ${map} > map.polisher.bam
      samtools sort map.polisher.bam -o map.polisher.sorted.bam
      samtools index map.polisher.sorted.bam
      """
    }
  }

/*
*------------------------------------------------------
*					POLISHING WITH RACON
*------------------------------------------------------
* If lrPolish = racon, start polishing the assembly using racon and long reads *
*/
	if (params.lrPolish == "racon" ){
		process raconLR {
			label 'raconLR'

			input:
			file reads from LongReadsPolisher_ch
			file assembly from AssemblyPolisher_ch
			file map from MapMinimap_ch

			output:
			file "assembly.racon${name}.fa" into PolisherAssembly_ch, AssemblyBuscoLR_ch
			file "${final_name}" optional true into FinalPolisherAssembly_ch

			script:
			iteration_polisherLR << "racon"
			name = iteration_polisherLR.size()
			final_name="assembly.racon"+LR_number+".fa"
			"""
			module load bioinfo/racon-v1.3.1
			racon -t 32 ${reads} ${map} ${assembly} > assembly.racon${name}.fa
			"""
		}
	}

/*
*------------------------------------------------------
*					POLISHING WITH WTDBG2
*------------------------------------------------------
* If LRPolish = wtdbg2, start polishing the assembly using wtdbg2 and long reads *
*/
	if (params.lrPolish == "wtdbg2"){
		process wtdbg2LR {
			label 'wtdbg2LR'

			input:
			file reads from LongReadsPolisher_ch
			file assembly from AssemblyPolisher_ch
			set map, index from MapSamtoolsLR_ch

			output:
			file "assembly.wtdbg2${name}.fa" into PolisherAssembly_ch, AssemblyBuscoLR_ch
			file "${final_name}" optional true into FinalPolisherAssembly_ch

			script:
			iteration_polisherLR << "wtdbg2"
			name = iteration_polisherLR.size()
			final_name="assembly.wtdbg2"+LR_number+".fa"
			"""
			module load bioinfo/wtdbg2-2.3
			samtools view -F0x900 ${map} | wtpoa-cns -t 32 -d ${assembly} -i - -fo assembly.wtdbg2${name}.fa
			"""
		}
	}

/*
*------------------------------------------------------
*					POLISHING WITH PILON
*------------------------------------------------------
* If LRPolish = pilon, start polishing the assembly using pilon and long reads *
*/
	if (params.lrPolish == "pilon"){
		process pilonLR {
			label 'pilonLR'

			input:
			file assembly from AssemblyPolisher_ch
			set map, index from MapSamtoolsLR_ch

			output:
			file "assembly.pilon${name}.fasta" into PolisherAssembly_ch, AssemblyBuscoLR_ch
			file "${final_name}" optional true into FinalPolisherAssembly_ch

			script:
			iteration_polisherLR << "pilon"
			name = iteration_polisherLR.size()
			final_name="assembly.pilon"+LR_number+".fasta"
			"""
			module load bioinfo/pilon-v1.22
			module load system/Java8
			java -Xmx32G -jar /usr/local/bioinfo/src/Pilon/pilon-v1.22/pilon-1.22.jar --genome ${assembly} --bam ${map} --fix bases,gaps --changes --output assembly.pilon${name}
			"""
		}
	}

/*
*------------------------------------------------------
*					POLISHING WITH MEDAKA
*------------------------------------------------------
* If LRPolish = medaka, start polishing the assembly using medaka and long reads *
*/
	if (params.lrPolish == "medaka"){
		process medakaLR {
			label 'medakaLR'

			input:
			file assembly from AssemblyPolisher_ch
			file reads from LongReadsPolisher_ch

			output:
			file "assembly.medaka${name}_output/consensus.fasta" into PolisherAssembly_ch, AssemblyBuscoLR_ch
			file "${final_name}" optional true into FinalPolisherAssembly_ch

			script:
			final_name="assembly.medaka"+LR_number+"_output/consensus.fasta"
			iteration_polisherLR << "medaka"
			name = iteration_polisherLR.size()
			"""
			export PS=\${PS:-''}
			export PS1=\${PS1:-''}
			source /usr/local/bioinfo/src/Medaka/medaka-0.5.0/venv/bin/activate
			module load bioinfo/medaka-0.5.0
			medaka_consensus -i ${reads} -d ${assembly} -o assembly.medaka${name}_output -t 4
			"""
		}
	}
}






/*
*=======================================================
* 				SHORT READS POLISHING PART
*=======================================================
//**********************
//* CONFIGURATION PART *
//**********************
* If Short reads are given *
*/
if (mode2) {
//* Create iteration condition, PolisherAssemblySR channel, split reads into multiple channels and the loop for multiple polish using short reads *
	ShortReads_ch.collect().into{ShortReadsAligner_ch; ShortReadsAligner2_ch; ShortReadsAligner3_ch; ShortReadsPolisher_ch; ShortReadsKat_ch}
  ShortReads_ch2.collect().into{ShortReadsAligner_ch2; ShortReadsAligner2_ch2; ShortReadsAligner3_ch2; ShortReadsPolisher_ch2; ShortReadsKat_ch2}
	iteration_polisherSR=[]
	PolisherAssemblySR_ch = Channel.create()

//* If no long reads detected, start first polishing with short reads, using initial given assembly *
	if (!mode1) {
		FinalPolisherAssembly_ch=Assembly_ch
	}

//* Create different condition if Params.NoChanges specified, to loop over process until the assembly don't change *
	if (params.noChanges){
		FinalPolisherAssembly_ch.mix(PolisherAssemblySR_ch)
							.into{AssemblyBwa_ch; AssemblyPolisher_ch}
		PolisherChangesSR_ch=Channel.fromPath('main.nf')
		PolisherChangesSR2_ch=Channel.create()
		PolisherChangesSR_ch.mix(PolisherChangesSR2_ch)
							.until{it.size()==0}
							.set{PolisherChangesSRFinal_ch}
	}
//* If we polish without using pilon *
  else if (params.srPolish != "pilon") {
		condition = { iteration_polisherSR.size()<SR_number ? it : Channel.STOP }
		FinalPolisherAssembly_ch.mix(PolisherAssemblySR_ch.map(condition))
							.into{AssemblyBwa_ch; AssemblyPolisher_ch; AssemblyPolisherVCF_ch; Assembly4Longranger_ch}
	}
//* If we polish with pilon without noChanges option *
  else {
		FinalPolisherAssembly_ch.into{AssemblyBwa_ch; AssemblyPolisher_ch; AssemblyDivide_ch; Assembly4Longranger_ch}
	}


//****************
//* PROCESS PART *
//****************
//* Map short reads to assembly using bwa-mem *
	if (params.shortReadsBwa1 && params.srPolish == "racon"){
    if (params.shortReadsBwa2){
		  process Dbwa_mem_sam {
			  label 'bwa_mem'

			  input:
			  file assembly from AssemblyBwa_ch
			  file reads1 from ShortReadsAligner_ch
        file reads2 from ShortReadsAligner_ch2

			  output:
			  file 'map.polisher.sam' into MapBwa_ch

			  script:
			  """
			  module load bioinfo/bwa-0.7.17
			  bwa index ${assembly}
			  bwa mem -t 8 ${assembly} ${reads1} ${reads2} >  map.polisher.sam
			  """
		  }
    }
    else {
      process bwa_mem_sam {
			  label 'bwa_mem'

			  input:
			  file assembly from AssemblyBwa_ch
			  file reads1 from ShortReadsAligner_ch

			  output:
			  file 'map.polisher.sam' into MapBwa_ch

			  script:
			  """
			  module load bioinfo/bwa-0.7.17
			  bwa index ${assembly}
			  bwa mem -t 8 ${assembly} ${reads1}  >  map.polisher.sam
			  """
		  }
    }
  }

  if (params.shortReadsBwa1 && params.srPolish != "racon"){
    if (params.shortReadsBwa2){
      process Dbwa_mem_bam {
        label 'bwa_samtools'
        input:
        file assembly from AssemblyBwa_ch
    		file reads1 from ShortReadsAligner_ch
        file reads2 from ShortReadsAligner_ch2

        output:
        set file("map.polisher.sorted.bam"), file("map.polisher.sorted.bam.bai") into MapSamtoolsSR_ch

        script:
        """
        module load bioinfo/bwa-0.7.17
        module load bioinfo/samtools-1.4
  			bwa index ${assembly}
  			bwa mem -t 8 ${assembly} ${reads1} ${reads2} | samtools view -S -b > map.polisher.bam
        samtools sort map.polisher.bam -o map.polisher.sorted.bam
        samtools index map.polisher.sorted.bam
        """
      }
    }
    else {
      process bwa_mem_bam {
        label 'bwa_samtools'
        input:
        file assembly from AssemblyBwa_ch
    		file reads1 from ShortReadsAligner_ch

        output:
        set file("map.polisher.sorted.bam"), file("map.polisher.sorted.bam.bai") into MapSamtoolsSR_ch

        script:
        """
        module load bioinfo/bwa-0.7.17
        module load bioinfo/samtools-1.4
  			bwa index ${assembly}
  			bwa mem -t 8 ${assembly} ${reads1} | samtools view -S -b > map.polisher.bam
        samtools sort map.polisher.bam -o map.polisher.sorted.bam
        samtools index map.polisher.sorted.bam
        """
      }
    }
  }

//* Map short reads to assembly using longranger *
  if (params.shortReadsLongranger && params.srPolish != "racon"){
    process mkref {
      label 'longranger_mkref'

      input:
      file assembly from Assembly4Longranger_ch

      output:
      file 'refdata*' into LongrangerRef_ch

      script:
      """
      module load bioinfo/longranger-2.2.2
      longranger mkref ${assembly}
      """
    }

    process align {
      label 'longranger_align'

      input:
      file reference from LongrangerRef_ch
      file reads from ShortReadsAligner_ch
      val samples from Sample

      output:
      set file ("ALIGN/outs/possorted_bam.bam"), file ("ALIGN/outs/possorted_bam.bam.bai") into MapSamtoolsSR_ch

      script:
      """
      module load bioinfo/longranger-2.2.2
      longranger align --reference=${reference} --id=ALIGN --fastqs=${reads} --localcores=20 --localmem=64 ${samples}
      """
    }
  }

/*
*------------------------------------------------------
*					POLISHING WITH PILON
*------------------------------------------------------
***************************************************************
//* If NoChanges option = false (Polishing with user options) *
***************************************************************
*/
	if (!params.noChanges){
//*****************************************************************************************
//* If SRPolish = pilon, start polishing the assembly using pilon and short reads         *
//* Maximum pilonSR polishing  = 3                                              				  *
//* If you want to polish more than 3 times with pilon, mail at clement.birbes@inra.fr    *
//* Or modify the code below by copy pasting LOOP THREE, add +1 to each numbered channels *
//* and create X more channel ligne 505                                                   *
//*****************************************************************************************

// LOOP ONE #################################################################################################################################
		if (params.srPolish == "pilon"){
			process divideContigs{
				label 'divideContigs'
				input:
				file assembly from AssemblyDivide_ch
				val length from chunckSize

				output:
				file "ctg.txt" into ListContigs_ch, ListContigs_ch2

				shell:
				'''
				module load system/Python-3.6.3
				FastaSplit.py --fasta !{assembly} --length !{length}
				'''
			}

      DividedContigs_ch = ListContigs_ch.splitCsv(sep:";")
			AssemblyBamIndex_ch = AssemblyPolisher_ch.combine(MapSamtoolsSR_ch)
			MergedInput_ch = DividedContigs_ch.combine(AssemblyBamIndex_ch)

			process pilonSR {
				label 'pilonSR'

				input:
				set num, contigs, assembly, bam, index from MergedInput_ch

				output:
				file "pilonSR*.fasta" into AssemblySplitSR_ch

				script:
				"""
				module load system/Java8
				java -Xmx32G -jar /usr/local/bioinfo/src/Pilon/pilon-v1.22/pilon-1.22.jar --genome ${assembly} --bam ${bam} --fix bases,gaps --changes --output pilonSR${num} --targets ${contigs}
				"""
			}

			process fasta_concat {
				label 'fasta_concat'

				input:
				file fastas from AssemblySplitSR_ch.collect()
        val length from chunckSize
        file ctg from ListContigs_ch2


				output:
				file "assembly.pilon${name}.fa" into PolisherAssemblySR2_ch, AssemblyBuscoSR1_ch
        file "${final_name}" optional true into AssemblySR1_ch

				script:
				iteration_polisherSR << "pilon"
				name=iteration_polisherSR.size()
				final_name = "assembly.pilon"+SR_number+".fa"
				"""
				module load system/Python-3.6.3
				FastaConcat.py --length ${length}
				mv pilonOut.fa assembly.pilon${name}.fa
				"""
			}

      if (params.allSteps){
        process buscoPilon1 {
    			label 'busco'

    			input:
    			file assembly from AssemblyBuscoSR1_ch
    			file lineage from LineageSR1_ch.collect()
    			val species from BUSCOspecies_ch

    			output:
    			file "run_BuscoSR${name}" into busco_first_assemblySR1_ch

    			script:
    			name = iteration_polisherSR.size()
    			"""
    			module load system/Python-3.6.3
    			module load bioinfo/augustus-3.3
    			module load bioinfo/busco-3.0.2
    			python3 /usr/local/bioinfo/src/BUSCO/busco-3.0.2/scripts/run_BUSCO.py -c 8 -i ${assembly} -l ${lineage} -m geno --limit 10 -o BuscoSR${name} -sp ${species}
    			"""
    		}
      }


//LOOP TWO #################################################################################################################################
			if (params.srNum >= 2) {
				PolisherAssemblySR2_ch.into{AssemblyBwa2_ch; Assembly4Longranger2_ch; AssemblyDivide2_ch; AssemblyPolisher2_ch}
				if (params.shortReadsBwa1){
          if (params.shortReadsBwa2){
            process Dbwa_mem_bam2 {
              label 'bwa_samtools'
              input:
              file assembly from AssemblyBwa2_ch
          		file reads1 from ShortReadsAligner2_ch
              file reads2 from ShortReadsAligner2_ch2

              output:
              set file("map.polisher.sorted.bam"), file("map.polisher.sorted.bam.bai") into MapPilon2_ch

              script:
              """
              module load bioinfo/bwa-0.7.17
              module load bioinfo/samtools-1.4
        			bwa index ${assembly}
        			bwa mem -t 8 ${assembly} ${reads1} ${reads2} | samtools view -S -b > map.polisher.bam
              samtools sort map.polisher.bam -o map.polisher.sorted.bam
              samtools index map.polisher.sorted.bam
              """
            }
          }
          else{
            process bwa_mem_bam2 {
              label 'bwa_samtools'
              input:
              file assembly from AssemblyBwa2_ch
          		file reads1 from ShortReadsAligner2_ch

              output:
              set file("map.polisher.sorted.bam"), file("map.polisher.sorted.bam.bai") into MapPilon2_ch

              script:
              """
              module load bioinfo/bwa-0.7.17
              module load bioinfo/samtools-1.4
        			bwa index ${assembly}
        			bwa mem -t 8 ${assembly} ${reads1} | samtools view -S -b > map.polisher.bam
              samtools sort map.polisher.bam -o map.polisher.sorted.bam
              samtools index map.polisher.sorted.bam
              """
            }
          }
        }

				if (params.shortReadsLongranger){
					process mkref2 {
						label 'longranger_mkref'

						input:
						file assembly from Assembly4Longranger2_ch

						output:
						file 'refdata*' into LongrangerRef2_ch

						script:
						"""
						module load bioinfo/longranger-2.2.2
						longranger mkref ${assembly}
						"""
					}

					process align2 {
						label 'longranger_align'

						input:
						file reference from LongrangerRef2_ch
						file reads from ShortReadsAligner2_ch
            val samples from Sample

						output:
						set file ("ALIGN/outs/possorted_bam.bam"), file ("ALIGN/outs/possorted_bam.bam.bai") into MapPilon2_ch

						script:
						"""
						module load bioinfo/longranger-2.2.2
						longranger align --maxjobs=16 --reference=${reference} --id=ALIGN --fastqs=${reads} --localcores=20 --localmem=64 ${samples}
						"""
					}
				}

				process divideContigs2{
					label 'divideContigs'
					input:
					file assembly from AssemblyDivide2_ch
					val length from chunckSize

					output:
					file "ctg.txt" into ListContigs2_ch, ListContigs2_ch2


					shell:
					'''
					module load system/Python-3.6.3
					FastaSplit.py --fasta !{assembly} --length !{length}
					'''
				}

        DividedContigs2_ch = ListContigs2_ch.splitCsv(sep:";")
  			AssemblyBamIndex2_ch = AssemblyPolisher2_ch.combine(MapPilon2_ch)
  			MergedInput2_ch = DividedContigs2_ch.combine(AssemblyBamIndex2_ch)

				process pilonSR2 {
					label 'pilonSR'

					input:
					set num, contigs, assembly, bam, index from MergedInput2_ch

					output:
					file "pilonSR*.fasta" into AssemblySplitSR2_ch

					script:
					"""
					module load system/Java8
					java -Xmx32G -jar /usr/local/bioinfo/src/Pilon/pilon-v1.22/pilon-1.22.jar --genome ${assembly} --bam ${bam} --fix bases,gaps --changes --output pilonSR${num} --targets ${contigs}
					"""
				}

				process fasta_concat2 {
					label 'fasta_concat'

					input:
					file fastas from AssemblySplitSR2_ch.collect()
          val length from chunckSize
          file ctg from ListContigs2_ch2


					output:
					file "assembly.pilon${name}.fa" into PolisherAssemblySR3_ch, AssemblyBuscoSR2_ch
          file "${final_name}" optional true into AssemblySR2_ch

					script:
					iteration_polisherSR << "pilon"
					name=iteration_polisherSR.size()
					final_name = "assembly.pilon"+SR_number+".fa"
					"""
					module load system/Python-3.6.3
					FastaConcat.py --length ${length}
					mv pilonOut.fa assembly.pilon${name}.fa
					"""
				}

        if (params.allSteps){
          process buscoPilon2 {
      			label 'busco'

      			input:
      			file assembly from AssemblyBuscoSR2_ch
      			file lineage from LineageSR2_ch.collect()
      			val species from BUSCOspecies_ch

      			output:
      			file "run_BuscoSR${name}" into busco_first_assemblySR2_ch

      			script:
      			name = iteration_polisherSR.size()
      			"""
      			module load system/Python-3.6.3
      			module load bioinfo/augustus-3.3
      			module load bioinfo/busco-3.0.2
      			python3 /usr/local/bioinfo/src/BUSCO/busco-3.0.2/scripts/run_BUSCO.py -c 8 -i ${assembly} -l ${lineage} -m geno --limit 10 -o BuscoSR${name} -sp ${species}
      			"""
      		}
        }


// LOOP THREE #################################################################################################################################
  			if (params.srNum >= 3) {
  				PolisherAssemblySR3_ch.into{AssemblyBwa3_ch; Assembly4Longranger3_ch; AssemblyDivide3_ch; AssemblyPolisher3_ch}
  				if (params.shortReadsBwa1){
            if (params.shortReadsBwa2){
              process Dbwa_mem_bam3 {
                label 'bwa_samtools'
                input:
                file assembly from AssemblyBwa3_ch
            		file reads1 from ShortReadsAligner3_ch
                file reads2 from ShortReadsAligner3_ch2

                output:
                set file("map.polisher.sorted.bam"), file("map.polisher.sorted.bam.bai") into MapPilon3_ch

                script:
                """
                module load bioinfo/bwa-0.7.17
                module load bioinfo/samtools-1.4
          			bwa index ${assembly}
          			bwa mem -t 8 ${assembly} ${reads1} ${reads2} | samtools view -S -b > map.polisher.bam
                samtools sort map.polisher.bam -o map.polisher.sorted.bam
                samtools index map.polisher.sorted.bam
                """
              }
            } else{
              process bwa_mem_bam3 {
                label 'bwa_samtools'
                input:
                file assembly from AssemblyBwa3_ch
            		file reads1 from ShortReadsAligner3_ch

                output:
                set file("map.polisher.sorted.bam"), file("map.polisher.sorted.bam.bai") into MapPilon3_ch

                script:
                """
                module load bioinfo/bwa-0.7.17
                module load bioinfo/samtools-1.4
          			bwa index ${assembly}
          			bwa mem -t 8 ${assembly} ${reads1} | samtools view -S -b > map.polisher.bam
                samtools sort map.polisher.bam -o map.polisher.sorted.bam
                samtools index map.polisher.sorted.bam
                """
              }
            }

          }

  				if (params.shortReadsLongranger){
  					process mkref3 {
  						label 'longranger_mkref'

  						input:
  						file assembly from Assembly4Longranger3_ch

  						output:
  						file 'refdata*' into LongrangerRef3_ch

  						script:
  						"""
  						module load bioinfo/longranger-2.2.2
  						longranger mkref ${assembly}
  						"""
  					}

  					process align3 {
  						label 'longranger_align'

  						input:
  						file reference from LongrangerRef3_ch
  						file reads from ShortReadsAligner3_ch
              val samples from Sample

  						output:
  						set file ("ALIGN/outs/possorted_bam.bam"), file ("ALIGN/outs/possorted_bam.bam.bai") into MapPilon3_ch

  						script:
  						"""
  						module load bioinfo/longranger-2.2.2
  						longranger align --maxjobs=16 --reference=${reference} --id=ALIGN --fastqs=${reads} --localcores=20 --localmem=64 ${samples}
  						"""
  					}
  				}

  				process divideContigs3{
  					label 'divideContigs'
  					input:
  					file assembly from AssemblyDivide3_ch
  					val length from chunckSize

  					output:
  					file "ctg.txt" into ListContigs3_ch, ListContigs3_ch2


  					shell:
  					'''
  					module load system/Python-3.6.3
  					FastaSplit.py --fasta !{assembly} --length !{length}
  					'''
  				}

          DividedContigs3_ch = ListContigs3_ch.splitCsv(sep:";")
    			AssemblyBamIndex3_ch = AssemblyPolisher3_ch.combine(MapPilon3_ch)
    			MergedInput3_ch = DividedContigs3_ch.combine(AssemblyBamIndex3_ch)

  				process pilonSR3 {
  					label 'pilonSR'

  					input:
  					set num, contigs, assembly, bam, index from MergedInput3_ch

  					output:
  					file "pilonSR*.fasta" into AssemblySplitSR3_ch

  					script:
  					"""
  					module load system/Java8
  					java -Xmx32G -jar /usr/local/bioinfo/src/Pilon/pilon-v1.22/pilon-1.22.jar --genome ${assembly} --bam ${bam} --fix bases,gaps --changes --output pilonSR${num} --targets ${contigs}
  					"""
  				}

  				process fasta_concat3 {
  					label 'fasta_concat'

  					input:
  					file fastas from AssemblySplitSR3_ch.collect()
            val length from chunckSize
            file ctg from ListContigs3_ch2


  					output:
  					file "assembly.pilon${name}.fa" into PolisherAssemblySR4_ch, AssemblyBuscoSR3_ch
            file "${final_name}" optional true into AssemblySR3_ch

  					script:
  					iteration_polisherSR << "pilon"
  					name=iteration_polisherSR.size()
  					final_name = "assembly.pilon"+SR_number+".fa"
  					"""
  					module load system/Python-3.6.3
  					FastaConcat.py --length ${length}
  					mv pilonOut.fa assembly.pilon${name}.fa
  					"""
  				}

          if (params.allSteps){
            process buscoPilon3 {
        			label 'busco'

        			input:
        			file assembly from AssemblyBuscoSR3_ch
        			file lineage from LineageSR3_ch.collect()
        			val species from BUSCOspecies_ch

        			output:
        			file "run_BuscoSR${name}" into busco_first_assemblySR3_ch

        			script:
        			name = iteration_polisherSR.size()
        			"""
        			module load system/Python-3.6.3
        			module load bioinfo/augustus-3.3
        			module load bioinfo/busco-3.0.2
        			python3 /usr/local/bioinfo/src/BUSCO/busco-3.0.2/scripts/run_BUSCO.py -c 8 -i ${assembly} -l ${lineage} -m geno --limit 10 -o BuscoSR${name} -sp ${species}
        			"""
        		}
          }
          AssemblySR3_ch.into{AssemblyKat_ch; AssemblyKat2_ch; AssemblyBuscoFinal_ch}
			  } else { AssemblySR2_ch.into{AssemblyKat_ch; AssemblyKat2_ch; AssemblyBuscoFinal_ch}}
      } else { AssemblySR1_ch.into{AssemblyKat_ch; AssemblyKat2_ch; AssemblyBuscoFinal_ch}}
    }


/*
*------------------------------------------------------
*					POLISHING WITH RACON
*------------------------------------------------------
* If SRPolish = racon, start polishing the assembly using racon and short reads *
*/
		if (params.srPolish == "racon"){
			process raconSR {
				label 'raconSR'

				input:
				file reads from ShortReadsPolisher_ch
				file assembly from AssemblyPolisher_ch
				file map from MapBwa_ch

				output:
				file "assembly.racon${name}.fa" into PolisherAssemblySR_ch,  AssemblyBuscoSR_ch
				file "${final_name}" optional true into AssemblyBuscoFinal_ch, AssemblyKat_ch, AssemblyKat2_ch

				script:
				iteration_polisherSR << "racon"
				name=iteration_polisherSR.size()
				final_name = "assembly.racon"+SR_number+".fa"
				"""
				module load bioinfo/racon-v1.3.1
				racon -t 32 ${reads} ${map} ${assembly} > assembly.racon${name}.fa
				"""
			}
		}

/*
*------------------------------------------------------
*					POLISHING WITH WTDBG2
*------------------------------------------------------
* If SRPolish = wtdbg2, start polishing the assembly using wtdbg2 and short reads *
*/
		if (params.srPolish == "wtdbg2"){
			process wtdbg2SR {
				label 'wtdbg2SR'

				input:
				file reads from ShortReadsPolisher_ch
				file assembly from AssemblyPolisher_ch
				set map, index from MapSamtoolsSR_ch

				output:
				file "assembly.wtdbg2${name}.fa" into PolisherAssemblySR_ch,  AssemblyBuscoSR_ch
				file "${final_name}" optional true into AssemblyBuscoFinal_ch, AssemblyKat_ch, AssemblyKat2_ch

				script:
				iteration_polisherSR << "wtdbg2"
				name=iteration_polisherSR.size()
				final_name = "assembly.wtdbg2"+SR_number+".fa"
				"""
				module load bioinfo/wtdbg2-2.3
				samtools view -F0x900 ${map} | wtpoa-cns -t 32 -d ${assembly} -x sam-sr -fo assembly.wtdbg2${name}.fa
				"""
			}
		}



/*
*------------------------------------------------------
*			POLISHING WITH FREEBAYES VCFTOOLS
*------------------------------------------------------
* If SRPolish = freebayes, start polishing the assembly using freebayes vcftools and short reads *
*/
		if (params.srPolish == "freebayes"){

			process freebayes {
				label 'freebayes'

				input:
				file assembly from AssemblyPolisher_ch
				set map, index from MapSamtoolsSR_ch

				output:
				file "assemblage.vcf" into vcf_ch

				script:
				"""
				module load bioinfo/freebayes-v1.2.0-2
				freebayes -f ${assembly} ${map} > assemblage.vcf
				"""
			}

      if (params.pattern != false){
			  process vcftoolsPattern {
				  label 'vcftools'

				  input:
				  file assembly from AssemblyPolisherVCF_ch
				  file vcf from vcf_ch

				  output:
				  file "assembly.freebayes${name}.fa" into PolisherAssemblySR_ch,  AssemblyBuscoSR_ch
          file "${final_name}" optional true into AssemblyBuscoFinal_ch, AssemblyKat_ch, AssemblyKat2_ch

				  script:
          iteration_polisherSR << "freebayes"
			  	name=iteration_polisherSR.size()
				  final_name = "assembly.freebayes"+SR_number+".fa"
			  	"""
				  module load bioinfo/tabix-0.2.5
			  	module load bioinfo/vcftools-0.1.15

				  grep -v ";AC=0;" ${vcf} | grep -v ";AN=0;" | grep -v ";AC=${stringPattern}" | grep -v ";AN=${stringPattern}" > assemblage_filtered.vcf
			  	bgzip assemblage_filtered.vcf
			  	tabix -p vcf assemblage_filtered.vcf.gz
				  cat ${assembly} | vcf-consensus assemblage_filtered.vcf.gz > assembly.freebayes${name}.fa
				  """
			  }
      }
      else {
        process vcftools {
				  label 'vcftools'

				  input:
				  file assembly from AssemblyPolisherVCF_ch
				  file vcf from vcf_ch

				  output:
				  file "assembly.freebayes${name}.fa" into PolisherAssemblySR_ch,  AssemblyBuscoSR_ch
          file "${final_name}" optional true into AssemblyBuscoFinal_ch, AssemblyKat_ch, AssemblyKat2_ch

				  script:
          iteration_polisherSR << "freebayes"
			  	name=iteration_polisherSR.size()
				  final_name = "assembly.freebayes"+SR_number+".fa"
			  	"""
				  module load bioinfo/tabix-0.2.5
			  	module load bioinfo/vcftools-0.1.15

				  grep -v ";AC=0;" ${vcf} | grep -v ";AN=0;" > assemblage_filtered.vcf
			  	bgzip assemblage_filtered.vcf
			  	tabix -p vcf assemblage_filtered.vcf.gz
				  cat ${assembly} | vcf-consensus assemblage_filtered.vcf.gz > assembly.freebayes${name}.fa
				  """
			  }
      }
		}
	}


/*
*------------------------------------------------------
*					POLISHING NoChanges
*------------------------------------------------------
* If NoChanges option = true (Polishing with pilon until no modification are declared) *
*/
	else {
		process pilonNoChanges {
			label 'pilonSR'

			input:
			file assembly from AssemblyPolisher_ch
			set map, index from MapSamtoolsSR_ch
			file change from PolisherChangesSRFinal_ch

			output:
			file "assembly.pilonSR${name}.fasta" into PolisherAssemblySR_ch, AssemblyBuscoSR_ch
			file "assembly.pilonSR${name}.changes" into PolisherChangesSR2_ch
			file "${final_name}" optional true into AssemblyBuscoFinal_ch, AssemblyKat_ch, AssemblyKat2_ch

			script:
			iteration_polisherSR << "pilon"
			name = iteration_polisherSR.size()
			final_name = "assembly.pilonSR"+SR_number+".fa"
			"""
			module load bioinfo/pilon-v1.22
			module load system/Java8

			java -Xmx32G -jar /usr/local/bioinfo/src/Pilon/pilon-v1.22/pilon-1.22.jar --genome ${assembly} --bam ${map} --fix bases,gaps --changes --output assembly.pilonSR${name}
			"""
		}
	}
}
/*
*=======================================================
* 					OUTPUT EVALUATION
*=======================================================
*------------------------------------------------------
*					BUSCO QUALITY
*------------------------------------------------------
* if lineage is specified, run busco quality after the last short reads polishing step
*/
if (mode3){
  if (!params.allSteps){
		process busco {
			label 'busco'

			input:
			file assembly from AssemblyBuscoFinal_ch
			file lineage from LineageSR_ch.collect()
			val species from BUSCOspecies_ch

			output:
			file "run_BuscoSR${name}" into busco_first_assemblyFinal_ch

			script:
			name = iteration_polisherSR.size()
			"""
			module load system/Python-3.6.3
			module load bioinfo/augustus-3.3
			module load bioinfo/busco-3.0.2
			python3 /usr/local/bioinfo/src/BUSCO/busco-3.0.2/scripts/run_BUSCO.py -c 8 -i ${assembly} -l ${lineage} -m geno --limit 10 -o BuscoSR${name} -sp ${species}
			"""
		}
  }

//* if params.allSteps = true run busco after each polishing steps
  if (params.allSteps && mode1){
		process buscoLRAllSteps {
			label 'busco'

			input:
			file assembly from AssemblyBuscoLR_ch
			file lineage from LineageLR_ch.collect()
			val species from BUSCOspecies_ch

			output:
			file "run_BuscoLR${name}" into busco_first_assemblyLR_ch

			script:
			name = iteration_polisherLR.size()
			"""
			module load system/Python-3.6.3
			module load bioinfo/augustus-3.3
			module load bioinfo/busco-3.0.2
			python3 /usr/local/bioinfo/src/BUSCO/busco-3.0.2/scripts/run_BUSCO.py -c 8 -i ${assembly} -l ${lineage} -m geno --limit 10 -o BuscoLR${name} -sp ${species}
			"""
		}
  }

  if (params.allSteps && mode2){
    AssemblyBuscoSR_ch=Channel.create()
    process buscoSRAllSteps {
			label 'busco'

			input:
			file assembly from AssemblyBuscoSR_ch
			file lineage from LineageSR_ch.collect()
			val species from BUSCOspecies_ch

			output:
			file "run_BuscoSR${name}" into busco_first_assemblySR_ch

			script:
			name = iteration_polisherSR.size()
			"""
			module load system/Python-3.6.3
			module load bioinfo/augustus-3.3
			module load bioinfo/busco-3.0.2
			python3 /usr/local/bioinfo/src/BUSCO/busco-3.0.2/scripts/run_BUSCO.py -c 8 -i ${assembly} -l ${lineage} -m geno --limit 10 -o BuscoSR${name} -sp ${species}
			"""
		}
  AssemblyBuscoSR_ch.close()
  }
}


//------------------------------------------------------
//				KAT QUALITY
//------------------------------------------------------
//* if params.kat=true run kat at the end of the pipeline (SR and/or LR vs Last available assembly generated by the pipeline.)
if (params.kat && mode2){
  if (params.shortReadsBwa2){
    process DkatSR {
      label 'katSR'

      input:
      file assembly from AssemblyKat_ch
      file reads1 from ShortReadsKat_ch
      file reads2 from ShortReadsKat_ch2

      output:
      file "*" into KatSROutput_ch

      script:
      """
      module load system/Miniconda3-4.4.10
      module load bioinfo/kat-2.4.1
      kat comp -t 8 '${reads1} ${reads2}' ${assembly}
      """
    }
  }
  else{
    process katSR {
      label 'katSR'

      input:
      file assembly from AssemblyKat_ch
      file reads1 from ShortReadsKat_ch

      output:
      file "*" into KatSROutput_ch

      script:
      """
      module load system/Miniconda3-4.4.10
      module load bioinfo/kat-2.4.1
      kat comp -t 8 ${reads1} ${assembly}
      """
    }
  }
}

if (params.kat && mode1){
  if (!mode2) {
		FinalPolisherAssembly_ch.set{AssemblyKat2_ch}
	}
  process katLR {
    label 'katLR'

    input:
    file assembly from AssemblyKat2_ch
    file reads from LongReadsKat_ch

    output:
    file "*" into KatLROutput_ch

    script:
    """
    module load system/Miniconda3-4.4.10
    module load bioinfo/kat-2.4.1
    kat comp -t 8 ${reads} ${assembly} -m 21 -p png
    """
  }
}
