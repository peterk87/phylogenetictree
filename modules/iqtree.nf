
process buildPhylogeneticIQTREE {
    
    tag { params.prefix }

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}"

    input:
    file(msa_mafft)

    output:
    path "*.iqtree", emit: ch_iqtree
    path "*.treefile", emit: ch_iqtree_newick

    script:
    """
    iqtree -s ${msa_mafft} -m GTR
    """
}


process visualizePhylogeneticGGTREE {

    tag { params.prefix }

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}"


    script:
    """
    phylogenetic_ggtree.r
    """
}