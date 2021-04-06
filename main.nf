#!/usr/bin/env nextflow
/*
========================================================================================
                         nf-core/phylogenetictree
========================================================================================
 nf-core/phylogenetictree Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/nf-core/phylogenetictree
----------------------------------------------------------------------------------------
*/

// enable dsl2
nextflow.preview.dsl = 2

include {catConsensusSequences; msaMAFFT} from './modules/msa'
include {buildPhylogeneticIQTREE; visualizePhylogeneticGGTREE; makeAlleles; rerootPhylogeneticTree; assignLineages} from './modules/iqtree'

// import subworkflows
include {articNcovNanopore} from './articNcovNanopore.nf' 

workflow {

    if (params.medaka|| params.nanopolish){
       // Check to see if we have barcodes
       nanoporeBarcodeDirs = file("${params.basecalled_fastq}/barcode*", type: 'dir', maxdepth: 1 )
       nanoporeNoBarcode = file("${params.basecalled_fastq}/*.fastq", type: 'file', maxdepth: 1) 

       if (nanoporeBarcodeDirs) {
            // Yes, barcodes!
            Channel.fromPath(nanoporeBarcodeDirs)
                .filter( ~/.*barcode[0-9]{1,4}$/ )
                .filter{ d ->
                            def count = 0
                            for (x in d.listFiles()) {
                                if (x.isFile()) {
                                    count += x.countFastq()
                                }
                            }
                            count > params.minReadsPerBarcode
                    }.set{ ch_fastqDirs }
                   
       }
       else if (nanoporeNoBarcode){
            // No, no barcodes
            Channel.fromPath( "${params.basecalled_fastq}", type: 'dir', maxDepth: 1 )
                    .set{ ch_fastqDirs }
       }
       else {
            println("Couldn't detect whether your Nanopore run was barcoded or not. Use --basecalled_fastq to point to the unmodified guppy output directory.")
            System.exit(1)
       }

        articNcovNanopore(ch_fastqDirs) 
    
    }


}