// ARTIC ncov workflow

// enable dsl2
nextflow.preview.dsl = 2

include {buildPhylogeneticIQTREE; visualizePhylogeneticGGTREE; makeAlleles; rerootPhylogeneticTree; assignLineages} from './modules/iqtree'

workflow phylogeneticTree {

    take:
      ch_runFastqDirs

    main:

      /*
      Multiuple Sequence Alignment
      */
      catConsensusSequences(articNcovNanopore.out.collect().combine(articDownloadScheme.out.reffasta))

      msaMAFFT(catConsensusSequences.out.ch_catconsensus)

      /*
      Build Phylogenetic Tree using IQTREE
      */

      buildPhylogeneticIQTREE(msaMAFFT.out.ch_msa_mafft)

      rerootPhylogeneticTree(buildPhylogeneticIQTREE.out.ch_iqtree_newick, articDownloadScheme.out.reffasta)

      makeAlleles(msaMAFFT.out.ch_msa_mafft)

      assignLineages(catConsensusSequences.out.ch_catconsensus)
    
      /*
      Visualization using bioconductor-ggtree 
      */

      visualizePhylogeneticGGTREE(rerootPhylogeneticTree.out.ch_reroot_iqtree, assignLineages.out.ch_lineage_report, makeAlleles.out.ch_alleles)

}


