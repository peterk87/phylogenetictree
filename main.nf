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


include {articDownloadScheme; articGuppyPlex; articMinIONMedaka; articRemoveUnmappedReads} from './modules/artic'
include {makeQCCSV; writeQCSummaryCSV} from './modules/qc'
include {collateSamples} from './modules/upload'
include {catConsensusSequences; msaMAFFT} from './modules/msa'
include {buildPhylogeneticIQTREE; visualizePhylogeneticGGTREE} from './modules/iqtree'



workflow {

    if (params.nanopore){
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
    // 
    /*
    STEP1: Run Artic workflow for nanopore sequencing
    */
    articDownloadScheme()
    
    articGuppyPlex(ch_fastqDirs.flatten())
    articMinIONMedaka(articGuppyPlex.out.fastq
                                      .filter{ it.countFastq() > params.minReadsArticGuppyPlex }
                                      .combine(articDownloadScheme.out.scheme))

    articRemoveUnmappedReads(articMinIONMedaka.out.mapped)

    makeQCCSV(articMinIONMedaka.out.ptrim.join(articMinIONMedaka.out.consensus_fasta, by: 0)
                           .combine(articDownloadScheme.out.reffasta))

    makeQCCSV.out.csv.splitCsv()
                       .unique()
                       .branch {
                           header: it[-1] == 'qc_pass'
                           fail: it[-1] == 'FALSE'
                           pass: it[-1] == 'TRUE'
                       }
                       .set { qc }

    writeQCSummaryCSV(qc.header.concat(qc.pass).concat(qc.fail).toList())
    
    collateSamples(qc.pass.map{ it[0] }
                           .join(articMinIONMedaka.out.consensus_fasta, by: 0)
                           .join(articRemoveUnmappedReads.out))
    

    /*
    STEP2: Multiple Sequence Alignment using MAFFT
    */

    catConsensusSequences(collateSamples.out.ch_samples_consensus.collect())

    msaMAFFT(catConsensusSequences.out.ch_catconsensus)

    /*
    STEP3: Build Phylogenetic Tree using IQTREE
    */

    buildPhylogeneticIQTREE(msaMAFFT.out.ch_msa_mafft)
    
    /*
    STEP4: Visualization using bioconductor-ggtree 
    */

    visualizePhylogeneticGGTREE()

    }


}