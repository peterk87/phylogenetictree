include {MSA_MAFFT} from './modules/msa'
include {BUILDPHYLOGENETIC_IQTREE; ASSIGNLINEAGES_PANGOLIN; REROOT_PHYLOGENETICTREE; DETERMINE_SNPS; PHYLOGENTICTREE_SNPS; PHYLOGENTICTREE_SHIPTV} from './modules/phylogenetictree'


workflow PHYLOGENETIC {
    take:
        ch_consensus_seqs
        ch_reference_fasta

    main:
        
        MSA_MAFFT(ch_consensus_seqs, ch_reference_fasta)

        BUILDPHYLOGENETIC_IQTREE(MSA_MAFFT.out.ch_msa_mafft)

        REROOT_PHYLOGENETICTREE(BUILDPHYLOGENETIC_IQTREE.out.ch_iqtree_newick, ch_reference_fasta)

        DETERMINE_SNPS(MSA_MAFFT.out.ch_msa_mafft)

        ASSIGNLINEAGES_PANGOLIN(ch_consensus_seqs)

        PHYLOGENTICTREE_SNPS(REROOT_PHYLOGENETICTREE.out.ch_reroot_iqtree, ASSIGNLINEAGES_PANGOLIN.out.ch_lineage_report, DETERMINE_SNPS.out.ch_alleles)

        PHYLOGENTICTREE_SHIPTV(REROOT_PHYLOGENETICTREE.out.ch_reroot_iqtree)
}