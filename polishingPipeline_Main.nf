#!/usr/bin/env nextflow
/*
========================================================================================
    			A S S E M B L Y   P O L I S H I N G   P I P E L I N E
========================================================================================
 Assembly Polishing Pipeline. Started October 2018.
 #### Authors ####
 Margot ZAHM <margot.zahm@inra.fr>
 Clément BIRBES <clement.birbes@inra.fr>
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
	    Assembly polishing Pipeline v0.1
	=========================================
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

    """.stripIndent()
}
//--reference			  Reference genome used for Quast comparison
//--genes           Gene and operon annotations used for Quast
//--poolseqSize     Size of the poolSeq sample
//--pattern         Maximum 0 for freebayes poolseq pattern (Eg: 3 = 0/0/0/1/..../1)

/*
*========================================================
*= 					SET UP CONFIGURATIONS				=
*========================================================
**********************************************
//* Show help message if --help is specified *
**********************************************
*/
params.help = false
if (params.help){
    helpMessage()
    exit 0
}

//*********************************
//* Default values initialization *
//*********************************
params.longReads		= false
params.shortReads		= false
params.assembly			= false

params.LRPolish 		= 'racon'
params.LRNum			  = 2
LR_number				    = params.LRNum
params.AlignerExt   = 'sam'

params.SRPolish 		= 'pilon'
params.SRNum		    = 2
SR_number				    = params.SRNum
params.NoChanges		= false

params.Outdir			  = './results/'

params.Lineage			= false
params.Species			= 'generic'

params.Kat          = false

//params.reference	= false
//params.genes			= false

params.Chunck			  = 10000000
chunckSize          = Channel.value(params.Chunck)

params.SRAligner		= 'longranger'

params.poolseqSize  = false
PSS                 = params.poolseqSize
params.pattern      = false
Patt                = params.pattern


//******************************************
//* Check if reads are given and not empty *
//******************************************
if (!params.longReads && !params.shortReads) {
	exit 1, "You must specify at least one short or long read file using --shortReads or --longReads."
}

if (params.longReads) {
	LongReads_ch=Channel.fromPath(params.longReads)
	                 .ifEmpty {exit 2, "Long Reads file not found: ${params.longReads}"}
	mode1=true
} else { mode1=false }

if (params.shortReads) {
	ShortReads_ch=Channel.fromPath(params.shortReads, type: 'dir')
				 .ifEmpty {exit 2, "Short Reads file not found: ${params.shortReads}"}
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
if (!params.Lineage) {
	mode3 = false
} else {
	Lineage_ch=Channel.fromPath(params.Lineage, type: 'dir')
						.ifEmpty {exit 1, "Assembly file not found: ${params.Lineage}"}
	BUSCOspecies_ch = Channel.value(params.Species)
	mode3 = true
}


//******************************************************
//* Check consistency of Long Reads Polishing argument *
//******************************************************
if (params.LRPolish != 'wtdbg2' && params.LRPolish != 'racon' && params.LRPolish != 'pilon' && params.LRPolish != 'medaka'){
	exit 1, "You must specify an available short reads polisher (wtdbg2, racon, medaka or pilon)"
}

if (params.LRNum.getClass() == Integer) {
	if  (params.LRNum > 5) {
		println "WARN : you are running ${params.raconNum} racon polishing"
	} else if (params.LRNum <=0) {
		exit 1, "Option --LRNum must be higher than 0"
	}
} else if (params.LRNum.getClass() != Integer) {
	exit 3, "Option --LRNum must be an integer : ${params.LRNum}"
}


//****************************************
//* If NoChanges option is NOT specified *
//****************************************
if (!params.NoChanges){
//*************************************
//* Check consistency of SR Polishing *
//*************************************
	if (params.SRPolish != 'wtdbg2' && params.SRPolish != 'racon' && params.SRPolish != 'pilon' && params.SRPolish != 'freebayes'){
		if (!params.NoChanges){
			exit 1, "You must specify an available short reads polisher (wtdbg2, racon, freebayes or pilon)"
		}
	}

	if (params.SRNum.getClass() == Integer) {
		if  (params.SRNum > 5) {
			println "WARN : you are running ${params.SRNum} pilon polishing"
		} else if (params.SRNum <=0) {
			exit 1, "Option --SRNum must be higher than 0"
		}
	} else if (params.SRNum.getClass() != Integer) {
		exit 3, "Option --SRNum must be an integer : ${params.SRNum}"
	}
}

//**************************
//* Chunck size evaluation *
//**************************
if (params.Chunck < 1000000){
	exit 1, "Using small chunck is not recommended, try to run pipeline with a chunk size bigger than 1,000,000"
}
if (params.Chunck > 100000000){
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


//*************************************************************
//* If Quast imput are specified A REVOIR PB VARIABLE GLOBALE *
//*************************************************************
//if (params.reference && params.genes) {
//	Ref=Channel.fromPath(params.reference)
//	                 .ifEmpty {exit 2, "Reference genome file not found: ${params.reference}"}
//	GeneAnnot=Channel.fromPath(params.genes)
//				.ifEmpty {exit 2, "Gene and operon annotations file not found: ${params.genes}"}
//	mode4=true
//} else { mode4=false }



/*
*========================================================
* 				LONG READS POLISHING PART
*========================================================
*********************************
//* If Long reads are specified *
*********************************
*/
if (mode1) {
//* Rename long reads file and add it to channels *
	process rename_long_reads {
		publishDir "${params.Outdir}"

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


//* Create iteration condition, PolisherAssembly channel and the loop for multiple polish using long reads *
	iteration_polisher=[]
	condition = { iteration_polisher.size()<LR_number ? it : Channel.STOP }

	PolisherAssembly_ch = Channel.create()

	Assembly_ch.mix( PolisherAssembly_ch.map(condition) )
			.into{AssemblyMinimap_ch; AssemblyPolisher_ch}


	if (params.LRPolish != "racon" || params.AlignerExt == "sam" ){
//* Map long reads to assembly using Minimap2 for Pilon, Medaka and wtdbg2 only *
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

  if (params.LRPolish == "racon" && params.AlignerExt == "paf" ) {
  //* Map long reads to assembly using Minimap2 for Racon only *
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

/*
*------------------------------------------------------
*					POLISHING WITH RACON
*------------------------------------------------------
**********************************************************************************
//* If LRPolish = racon, start polishing the assembly using racon and long reads *
**********************************************************************************
*/
	if (params.LRPolish == "racon" ){
		process raconLR {
			label 'raconLR'

			input:
			file reads from LongReadsPolisher_ch
			file assembly from AssemblyPolisher_ch
			file map from MapMinimap_ch

			output:
			file "assembly.racon${name}.fa" into PolisherAssembly_ch
			file "${final_name}" optional true into FinalPolisherAssembly_ch

			script:
			iteration_polisher << "racon"
			name = iteration_polisher.size()
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
************************************************************************************
//* If LRPolish = wtdbg2, start polishing the assembly using wtdbg2 and long reads *
************************************************************************************
*/
	if (params.LRPolish == "wtdbg2"){
		process wtdbg2LR {
			label 'wtdbg2LR'

			input:
			file reads from LongReadsPolisher_ch
			file assembly from AssemblyPolisher_ch
			file map from MapMinimap_ch

			output:
			file "assembly.wtdbg2${name}.fa" into PolisherAssembly_ch
			file "${final_name}" optional true into FinalPolisherAssembly_ch

			script:
			iteration_polisher << "wtdbg2"
			name = iteration_polisher.size()
			final_name="assembly.wtdbg2"+LR_number+".fa"
			"""
			module load bioinfo/samtools-1.9
			module load bioinfo/wtdbg2-2.3

			samtools sort -@4 ${map} > sorted.bam
			samtools view -F0x900 sorted.bam | wtpoa-cns -t 32 -d ${assembly} -i - -fo assembly.wtdbg2${name}.fa
			"""
		}
	}

/*
*------------------------------------------------------
*					POLISHING WITH PILON
*------------------------------------------------------
**********************************************************************************
//* If LRPolish = pilon, start polishing the assembly using pilon and long reads *
**********************************************************************************
*/
	if (params.LRPolish == "pilon"){
		process pilonLR {
			label 'pilonLR'

			input:
			file assembly from AssemblyPolisher_ch
			file map from MapMinimap_ch

			output:
			file "assembly.pilon${name}.fa" into PolisherAssembly_ch
			file "${final_name}" optional true into FinalPolisherAssembly_ch

			script:
			iteration_polisher << "pilon"
			name = iteration_polisher.size()
			final_name="assembly.pilon"+LR_number+".fa"
			"""
			module load bioinfo/samtools-1.4
			module load bioinfo/pilon-v1.22
			module load system/Java8

			samtools view -S -b ${map} > map.polisher.bam
			samtools sort map.polisher.bam -o map.polisher.sorted.bam
			samtools index map.polisher.sorted.bam
			java -Xmx32G -jar /usr/local/bioinfo/src/Pilon/pilon-v1.22/pilon-1.22.jar --genome ${assembly} --bam map.polisher.sorted.bam --fix bases,gaps --changes --output pilon${name}
			"""
		}
	}

/*
*------------------------------------------------------
*					POLISHING WITH MEDAKA
*------------------------------------------------------
************************************************************************************
//* If LRPolish = medaka, start polishing the assembly using medaka and long reads *
************************************************************************************
*/
	if (params.LRPolish == "medaka"){
		process medakaLR {
			label 'medakaLR'

			input:
			file assembly from AssemblyPolisher_ch
			file reads from LongReadsPolisher_ch

			output:
			file "assembly.medaka${name}_output/consensus.fasta" into PolisherAssembly_ch
			file "${final_name}" optional true into FinalPolisherAssembly_ch

			script:
			final_name="assembly.medaka"+LR_number+"_output/consensus.fasta"
			iteration_polisher << "medaka"
			name = iteration_polisher.size()
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
*********************************
//* If Long reads are specified *
*********************************
*/
if (mode2) {
	ShortReads_ch.collect().into{ShortReadsAligner_ch; ShortReadsAligner2_ch; ShortReadsAligner3_ch; ShortReadsPolisher_ch; ShortReadsQuast_ch; ShortReadsKat_ch}

//* Create iteration condition, PolisherAssemblySR channel and the loop for multiple polish using short reads *
//* Different condition if Params.NoChanges specified, to loop over process until the assembly don't change *
	iteration_polisherSR=[]
	PolisherAssemblySR_ch = Channel.create()

//* If no long reads detected, starting first polishing with short reads *
	if (!mode1) {
		FinalPolisherAssembly_ch=Assembly_ch
	}

	if (params.NoChanges){
		FinalPolisherAssembly_ch.mix(PolisherAssemblySR_ch)
							.into{AssemblyBwa_ch; AssemblyPolisher_ch}
		PolisherChangesSR_ch=Channel.fromPath('Main.nf')
		PolisherChangesSR2_ch=Channel.create()
		PolisherChangesSR_ch.mix(PolisherChangesSR2_ch)
							.until{it.size()==0}
							.set{PolisherChangesSRFinal_ch}

	} else if (params.SRPolish != "pilon") {
		condition = { iteration_polisherSR.size()<SR_number ? it : Channel.STOP }
		FinalPolisherAssembly_ch.mix(PolisherAssemblySR_ch.map(condition))
							.into{AssemblyBwa_ch; AssemblyPolisher_ch; AssemblyPolisherBis_ch; AssemblyDivide_ch; Assembly4Longranger_ch}
	} else {
		FinalPolisherAssembly_ch.into{AssemblyBwa_ch; AssemblyPolisher_ch; AssemblyPolisherBis_ch; AssemblyDivide_ch; Assembly4Longranger_ch}
	}


//* Map short reads to assembly using bwa-mem *
	if (params.SRAligner == "samtools"){
		process bwa_mem {
			label 'bwa_mem'

			input:
			file assembly from AssemblyBwa_ch
			file reads from ShortReadsAligner_ch

			output:
			file 'map.polisher.sam' into MapBwa_ch

			script:
			"""
			module load bioinfo/samtools-1.4
			module load bioinfo/bwa-0.7.17
			bwa index ${assembly}
			bwa mem -t 8 ${assembly} ${reads} >  map.polisher.sam

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
	if (!params.NoChanges){
//*********************************************************************************
//* If SRPolish = pilon, start polishing the assembly using pilon and short reads *
//* One block for each polishing. Impossible to loop + parallelize				  *
//*********************************************************************************

// LOOP ONE #################################################################################################################################
		if (params.SRPolish == "pilon"){
			if (params.SRAligner == "samtools"){
				process samtools{
					label 'samtools'
					input:
					file map from MapBwa_ch

					output:
					set file("map.polisher.sorted.bam"), file("map.polisher.sorted.bam.bai") into MapPilon_ch

					script:
					"""
					module load bioinfo/samtools-1.4
					samtools view -S -b ${map} > map.polisher.bam
					samtools sort map.polisher.bam -o map.polisher.sorted.bam
					samtools index map.polisher.sorted.bam
					"""
				}
			}

			if (params.SRAligner == "longranger"){
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

					output:
					set file ("ALIGN/outs/possorted_bam.bam"), file ("ALIGN/outs/possorted_bam.bam.bai") into MapPilon_ch

					script:
					"""
					module load bioinfo/longranger-2.2.2
					longranger align --maxjobs=16 --reference=${reference} --id=ALIGN --fastqs=${reads}
					"""
				}
			}

			process divideContigs{
				label 'divideContigs'
				input:
				file assembly from AssemblyDivide_ch
				val length from chunckSize

				output:
				file "ctg.txt" into ListContigs_ch

				shell:
				'''
				module load system/Python-3.6.3
				python3.6 ../../../FastaSplit.py --fasta !{assembly} --length !{length}
				'''
			}

      DividedContigs_ch = ListContigs_ch.splitCsv(sep:";")
			AssemblyBamIndex_ch = AssemblyPolisher_ch.combine(MapPilon_ch)
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
				java -Xmx32G -jar /usr/local/bioinfo/src/Pilon/pilon-v1.22/pilon-1.22.jar --genome ${assembly} --bam ${bam} --fix bases,gaps --changes --output pilonSR${contigs} --targets ${contigs}
				"""
			}

			process fasta_concat {
				label 'fasta_concat'

				input:
				file fastas from AssemblySplitSR_ch.collect()
        val length from chunckSize

				output:
				file "assembly.pilon${name}.fa" into PolisherAssemblySR2_ch

				script:
				iteration_polisherSR << "pilon"
				name=iteration_polisherSR.size()
				final_name = "assembly.pilon"+SR_number+".fa"
				"""
				module load system/Python-3.6.3
				python3.6 ../../../FastaConcat.py --length ${length}
				mv pilonOut.fa assembly.pilon${name}.fa
				"""
			}


//LOOP TWO #################################################################################################################################
			if (params.SRNum >= 2) {
				PolisherAssemblySR2_ch.into{AssemblyBwa2_ch; Assembly4Longranger2_ch; AssemblyDivide2_ch; AssemblyPolisher2_ch}
				if (params.SRAligner == "samtools"){
					process bwa_mem2 {
						label 'bwa_mem'

						input:
						file assembly from AssemblyBwa2_ch
						file reads from ShortReadsAligner2_ch

						output:
						file 'map.polisher.sam' into MapBwa2_ch

						script:
						"""
						module load bioinfo/samtools-1.4
						module load bioinfo/bwa-0.7.17
						bwa index ${assembly}
						bwa mem -t 8 ${assembly} ${reads} >  map.polisher.sam

						"""
					}

					process samtools2{
						label 'samtools'
						input:
						file map from MapBwa2_ch

						output:
						set file("map.polisher.sorted.bam"), file("map.polisher.sorted.bam.bai") into MapPilon2_ch

						script:
						"""
						module load bioinfo/samtools-1.4
						samtools view -S -b ${map} > map.polisher.bam
						samtools sort map.polisher.bam -o map.polisher.sorted.bam
						samtools index map.polisher.sorted.bam
						"""
					}
        }

				if (params.SRAligner == "longranger"){
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

						output:
						set file ("ALIGN/outs/possorted_bam.bam"), file ("ALIGN/outs/possorted_bam.bam.bai") into MapPilon2_ch

						script:
						"""
						module load bioinfo/longranger-2.2.2
						longranger align --maxjobs=16 --reference=${reference} --id=ALIGN --fastqs=${reads}
						"""
					}
				}

				process divideContigs2{
					label 'divideContigs'
					input:
					file assembly from AssemblyDivide2_ch
					val length from chunckSize

					output:
					file "ctg.txt" into ListContigs2_ch

					shell:
					'''
					module load system/Python-3.6.3
					python3.6 ../../../FastaSplit.py --fasta !{assembly} --length !{length}
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
					java -Xmx32G -jar /usr/local/bioinfo/src/Pilon/pilon-v1.22/pilon-1.22.jar --genome ${assembly} --bam ${bam} --fix bases,gaps --changes --output pilonSR${contigs} --targets ${contigs}
					"""
				}

				process fasta_concat2 {
					label 'fasta_concat'

					input:
					file fastas from AssemblySplitSR2_ch.collect()
          val length from chunckSize

					output:
					file "assembly.pilon${name}.fa" into PolisherAssemblySR3_ch

					script:
					iteration_polisherSR << "pilon"
					name=iteration_polisherSR.size()
					final_name = "assembly.pilon"+SR_number+".fa"
					"""
					module load system/Python-3.6.3
					python3.6 ../../../FastaConcat.py --length ${length}
					mv pilonOut.fa assembly.pilon${name}.fa
					"""
				}


// LOOP THREE #################################################################################################################################
  			if (params.SRNum >= 3) {
  				PolisherAssemblySR3_ch.into{AssemblyBwa3_ch; Assembly4Longranger3_ch; AssemblyDivide3_ch; AssemblyPolisher3_ch}
  				if (params.SRAligner == "samtools"){
  					process bwa_mem3 {
  						label 'bwa_mem'

  						input:
  						file assembly from AssemblyBwa3_ch
  						file reads from ShortReadsAligner3_ch

  						output:
  						file 'map.polisher.sam' into MapBwa3_ch

  						script:
  						"""
  						module load bioinfo/samtools-1.4
  						module load bioinfo/bwa-0.7.17
  						bwa index ${assembly}
  						bwa mem -t 8 ${assembly} ${reads} >  map.polisher.sam

  						"""
  					}

  					process samtools3{
  						label 'samtools'
  						input:
  						file map from MapBwa3_ch

  						output:
  						set file("map.polisher.sorted.bam"), file("map.polisher.sorted.bam.bai") into MapPilon3_ch

  						script:
  						"""
  						module load bioinfo/samtools-1.4
  						samtools view -S -b ${map} > map.polisher.bam
  						samtools sort map.polisher.bam -o map.polisher.sorted.bam
  						samtools index map.polisher.sorted.bam
  						"""
  					}
          }

  				if (params.SRAligner == "longranger"){
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

  						output:
  						set file ("ALIGN/outs/possorted_bam.bam"), file ("ALIGN/outs/possorted_bam.bam.bai") into MapPilon3_ch

  						script:
  						"""
  						module load bioinfo/longranger-2.2.2
  						longranger align --maxjobs=16 --reference=${reference} --id=ALIGN --fastqs=${reads}
  						"""
  					}
  				}

  				process divideContigs3{
  					label 'divideContigs'
  					input:
  					file assembly from AssemblyDivide3_ch
  					val length from chunckSize

  					output:
  					file "ctg.txt" into ListContigs3_ch

  					shell:
  					'''
  					module load system/Python-3.6.3
  					python3.6 ../../../FastaSplit.py --fasta !{assembly} --length !{length}
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
  					java -Xmx32G -jar /usr/local/bioinfo/src/Pilon/pilon-v1.22/pilon-1.22.jar --genome ${assembly} --bam ${bam} --fix bases,gaps --changes --output pilonSR${contigs} --targets ${contigs}
  					"""
  				}

  				process fasta_concat3 {
  					label 'fasta_concat'

  					input:
  					file fastas from AssemblySplitSR3_ch.collect()
            val length from chunckSize

  					output:
  					file "assembly.pilon${name}.fa" into PolisherAssemblySR4_ch

  					script:
  					iteration_polisherSR << "pilon"
  					name=iteration_polisherSR.size()
  					final_name = "assembly.pilon"+SR_number+".fa"
  					"""
  					module load system/Python-3.6.3
  					python3.6 ../../../FastaConcat.py --length ${length}
  					mv pilonOut.fa assembly.pilon${name}.fa
  					"""
  				}
          PolisherAssemblySR4_ch.into{AssemblyBuscoSR_ch; AssemblyKat_ch; AssemblyKat2_ch}
			  } else { PolisherAssemblySR3_ch.into{AssemblyBuscoSR_ch; AssemblyKat_ch; AssemblyKat2_ch}}
      } else { PolisherAssemblySR2_ch.into{AssemblyBuscoSR_ch; AssemblyKat_ch; AssemblyKat2_ch}}
    }


//##############################################################################################################################
//# To get more polishing with pilon, copy one entire loop and change output and input channels (add +1 to each numerotations) #
//##############################################################################################################################

/*
*------------------------------------------------------
*					POLISHING WITH RACON
*------------------------------------------------------
***********************************************************************************
//* If SRPolish = racon, start polishing the assembly using racon and short reads *
***********************************************************************************
*/
		if (params.SRPolish == "racon"){
			process raconSR {
				label 'raconSR'

				input:
				file reads from ShortReadsPolisher_ch
				file assembly from AssemblyPolisher_ch
				file map from MapBwa_ch

				output:
				file "assembly.racon${name}.fa" into PolisherAssemblySR_ch
				file "${final_name}" optional true into AssemblyBuscoSR_ch, AssemblyKat_ch, AssemblyKat2_ch

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
*************************************************************************************
//* If SRPolish = wtdbg2, start polishing the assembly using wtdbg2 and short reads *
*************************************************************************************
*/
		if (params.SRPolish == "wtdbg2"){
			process wtdbg2SR {
				label 'wtdbg2SR'

				input:
				file reads from ShortReadsPolisher_ch
				file assembly from AssemblyPolisher_ch
				file map from MapBwa_ch

				output:
				file "assembly.wtdbg2${name}.fa" into PolisherAssemblySR_ch
				file "${final_name}" optional true into AssemblyBuscoSR_ch, AssemblyKat_ch, AssemblyKat2_ch

				script:
				iteration_polisherSR << "wtdbg2"
				name=iteration_polisherSR.size()
				final_name = "assembly.wtdbg2"+SR_number+".fa"
				"""
				module load bioinfo/samtools-1.9
				module load bioinfo/wtdbg2-2.3

				samtools sort -@4 ${map} > sorted.bam
				samtools view -F0x900 sorted.bam | wtpoa-cns -t 32 -d ${assembly} -x sam-sr -fo assembly.wtdbg2${name}.fa
				"""
			}
		}



/*
*------------------------------------------------------
*			POLISHING WITH FREEBAYES VCFTOOLS
*------------------------------------------------------
****************************************************************************************************
//* If SRPolish = freebayes, start polishing the assembly using freebayes vcftools and short reads *
****************************************************************************************************
*/
		if (params.SRPolish == "freebayes"){

			process samtoolsFreebayes{
				label 'samtools'
				input:
				file map from MapBwa_ch

				output:
				set file("map.polisher.sorted.bam"), file("map.polisher.sorted.bam.bai") into MapFreeBayes_ch

				script:
				"""
				module load bioinfo/samtools-1.4
				samtools view -S -b ${map} > map.polisher.bam
				samtools sort map.polisher.bam -o map.polisher.sorted.bam
				samtools index map.polisher.sorted.bam
				"""
			}

			process freebayes {
				label 'freebayes'

				input:
				file assembly from AssemblyPolisher_ch
				set map, index from MapFreeBayes_ch

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
				  file assembly from AssemblyPolisherBis_ch
				  file vcf from vcf_ch

				  output:
				  file "assembly.freebayes${name}.fa" into PolisherAssemblySR_ch
          file "${final_name}" optional true into AssemblyBuscoSR_ch, AssemblyKat_ch, AssemblyKat2_ch

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
      } else {
        process vcftools {
				  label 'vcftools'

				  input:
				  file assembly from AssemblyPolisherBis_ch
				  file vcf from vcf_ch

				  output:
				  file "assembly.freebayes${name}.fa" into PolisherAssemblySR_ch
          file "${final_name}" optional true into AssemblyBuscoSR_ch, AssemblyKat_ch, AssemblyKat2_ch

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
******************************************************************************************
//* If NoChanges option = true (Polishing with pilon until no modification are declared) *
******************************************************************************************
*/
	else {
		process pilonNoChanges {
			label 'pilonSR'

			input:
			file assembly from AssemblyPolisher_ch
			file map from MapBwa_ch
			file change from PolisherChangesSRFinal_ch

			output:
			file "pilonSR${name}.fasta" into PolisherAssemblySR_ch
			file "pilonSR${name}.changes" into PolisherChangesSR2_ch
			file "${final_name}" optional true into AssemblyBuscoSR_ch, AssemblyKat_ch, AssemblyKat2_ch

			script:
			iteration_polisherSR << "pilon"
			name = iteration_polisherSR.size()
			final_name = "pilonSR"+SR_number+".fasta"
			"""
			module load bioinfo/samtools-1.4
			module load bioinfo/pilon-v1.22
			module load system/Java8

			samtools view -S -b ${map} > map.polisher.bam
			samtools sort map.polisher.bam -o map.polisher.sorted.bam
			samtools index map.polisher.sorted.bam
			java -Xmx32G -jar /usr/local/bioinfo/src/Pilon/pilon-v1.22/pilon-1.22.jar --genome ${assembly} --bam map.polisher.sorted.bam --fix bases,gaps --changes --output pilonSR${name}
			"""
		}
	}

/*
*=======================================================
* 					OUTPUT EVALUATION
*=======================================================
*------------------------------------------------------
*					BUSCO QUALITY
*------------------------------------------------------
*/
	if (mode3){
		process buscoSR {
			label 'buscoSR'

			input:
			file assembly from AssemblyBuscoSR_ch
			file lineage from Lineage_ch.collect()
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
	}
//------------------------------------------------------
//				KAT QUALITY
//------------------------------------------------------

  if (params.Kat && mode2){
    process katSR {
      label 'katSR'

      input:
      file assembly from AssemblyKat_ch
      file reads from ShortReadsKat_ch

      output:
      file "*" into KatSROutput_ch

      script:
      """
      module load system/Miniconda3-4.4.10
      module load bioinfo/kat-2.4.1
      kat comp -t 8 ${reads} ${assembly}
      """
    }
  }

  if (params.Kat && mode1){
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
      kat comp -t 8 ${reads} ${assembly}
      """
    }
  }



/*
*------------------------------------------------------
* QUAST QUALITY (problème variables globale)
*------------------------------------------------------
*/
// Decommenter + rajouter output AssemblyQuast dans tous les process pour reutiliser
//		process quast {
//			label 'quast'
//
//			input:
//			file assembly from AssemblyQuast
//			file reference from Ref
//			file genes from GeneAnnot
//			file reads from ShortReadsQuast
//
//			output:
//			file "quast${name}_output" into Quast_ComparisonSR
//
//			script:
//			name = iteration_polisherSR.size()
//			"""
//			source /usr/local/bioinfo/src/QUAST/quast-5.0.0rc1_env/bin/activate
//			export _OLD_VIRTUAL_PATH=\${_OLD_VIRTUAL_PATH:-''}
//			export _OLD_VIRTUAL_PYTHONHOME=\${_OLD_VIRTUAL_PYTHONHOME:-''}
//			export _OLD_VIRTUAL_PS1=\${_OLD_VIRTUAL_PS1:-''}
//			export ZSH_VERSION=\${ZSH_VERSION:-''}
//			export VIRTUAL_ENV_DISABLE_PROMPT=\${VIRTUAL_ENV_DISABLE_PROMPT:-''}
//			module load system/R-3.4.3
//			module load bioinfo/quast-5.0.0rc1
//			quast.py ${assembly} -R ${reference} -g ${genes} --single ${reads} -o quast${name}_output
//			"""
//		}
//	}
}
