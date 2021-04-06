// ARTIC ncov workflow

// enable dsl2
nextflow.preview.dsl = 2

// import modules
include {articDownloadScheme} from './modules/artic.nf' 
include {articGuppyPlex} from './modules/artic.nf' 
include {articMinIONNanopolish} from  './modules/artic.nf' 
include {articMinIONMedaka} from  './modules/artic.nf'
include {articRemoveUnmappedReads} from './modules/artic.nf' 
include {makeQCCSV} from './modules/qc.nf'
include {writeQCSummaryCSV} from './modules/qc.nf'
include {collateSamples} from './modules/upload.nf'

include {catConsensusSequences; msaMAFFT} from './modules/msa'
include {buildPhylogeneticIQTREE; visualizePhylogeneticGGTREE; makeAlleles; rerootPhylogeneticTree; assignLineages} from './modules/iqtree'


// workflow component for artic pipeline
workflow sequenceAnalysisNanopolish {
    take:
      ch_runFastqDirs
      ch_fast5Pass
      ch_seqSummary
    
    main:
      articDownloadScheme()
      
      articGuppyPlex(ch_runFastqDirs.flatten())

      articMinIONNanopolish(articGuppyPlex.out.fastq
                                          .filter{ it.countFastq() > params.minReadsArticGuppyPlex }
                                          .combine(articDownloadScheme.out.scheme)
                                          .combine(ch_fast5Pass)
                                          .combine(ch_seqSummary))

      articRemoveUnmappedReads(articMinIONNanopolish.out.mapped)

      makeQCCSV(articMinIONNanopolish.out.ptrim
                                     .join(articMinIONNanopolish.out.consensus_fasta, by: 0)
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
                           .join(articMinIONNanopolish.out.consensus_fasta, by: 0)
                           .join(articRemoveUnmappedReads.out))
      /*
      Multiple Sequence Alignment
      */
      catConsensusSequences(collateSamples.out.ch_samples_consensus.collect().combine(articDownloadScheme.out.reffasta))
      msaMAFFT(catConsensusSequences.out.ch_catconsensus)

      /*
      Build Phylogenetic Tree
      */
      buildPhylogeneticIQTREE(msaMAFFT.out.ch_msa_mafft)
      rerootPhylogeneticTree(buildPhylogeneticIQTREE.out.ch_iqtree_newick, articDownloadScheme.out.reffasta)
      makeAlleles(msaMAFFT.out.ch_msa_mafft)
      assignLineages(catConsensusSequences.out.ch_catconsensus)
      
      /*
      Visualize Phylogenetic Tree
      */
      visualizePhylogeneticGGTREE(rerootPhylogeneticTree.out.ch_reroot_iqtree, assignLineages.out.ch_lineage_report, makeAlleles.out.ch_alleles)

    emit:
      qc_pass = collateSamples.out
      reffasta = articDownloadScheme.out.reffasta
      vcf = articMinIONNanopolish.out.vcf

}

workflow sequenceAnalysisMedaka {
    take:
      ch_runFastqDirs

    main:
      articDownloadScheme()

      articGuppyPlex(ch_runFastqDirs.flatten())

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
      Multiple Sequence Alignment
      */
      catConsensusSequences(collateSamples.out.ch_samples_consensus.collect().combine(articDownloadScheme.out.reffasta))
      msaMAFFT(catConsensusSequences.out.ch_catconsensus)

      /*
      Build Phylogenetic Tree
      */
      buildPhylogeneticIQTREE(msaMAFFT.out.ch_msa_mafft)
      rerootPhylogeneticTree(buildPhylogeneticIQTREE.out.ch_iqtree_newick, articDownloadScheme.out.reffasta)
      makeAlleles(msaMAFFT.out.ch_msa_mafft)
      assignLineages(catConsensusSequences.out.ch_catconsensus)

      /*
      Visualize Phylogenetic Tree
      */
      visualizePhylogeneticGGTREE(rerootPhylogeneticTree.out.ch_reroot_iqtree, assignLineages.out.ch_lineage_report, makeAlleles.out.ch_alleles)

    emit:
      qc_pass = collateSamples.out
      reffasta = articDownloadScheme.out.reffasta
      vcf = articMinIONMedaka.out.vcf

}


workflow articNcovNanopore {

    take:
      ch_fastqDirs    
      
    main:
    
      if ( params.nanopolish ) {
          Channel.fromPath( "${params.fast5_pass}" )
                 .set{ ch_fast5Pass }

          Channel.fromPath( "${params.sequencing_summary}" )
                 .set{ ch_seqSummary }

          sequenceAnalysisNanopolish(ch_fastqDirs, ch_fast5Pass, ch_seqSummary)

          sequenceAnalysisNanopolish.out.vcf.set{ ch_nanopore_vcf }

          sequenceAnalysisNanopolish.out.reffasta.set{ ch_nanopore_reffasta }
      
      } else if (params.medaka) {

          sequenceAnalysisMedaka(ch_fastqDirs)

          sequenceAnalysisMedaka.out.vcf.set{ ch_nanopore_vcf }

          sequenceAnalysisMedaka.out.reffasta.set{ ch_nanopore_reffasta }
          
      }
}   

