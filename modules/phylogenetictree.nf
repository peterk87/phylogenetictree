process BUILDPHYLOGENETIC_IQTREE {
    
    tag { params.prefix }

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}"

    input:
    file(msa_mafft)

    output:
    path "*.iqtree", emit: ch_iqtree
    path "*.treefile", emit: ch_iqtree_newick

    script:
    """
    iqtree -s ${msa_mafft} -m ${params.substitution_model}
    """
}


process ASSIGNLINEAGES {
    
    tag { params.prefix }

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}"

    input:
    file(cat_consensus_seqeuences)

    output:
    path "*.csv", emit: ch_lineage_report


    script:
    """
    pangolin -t 2 --outfile lineage_report.csv ${cat_consensus_seqeuences}
    """
}

process REROOT_PHYLOGENETICTREE {
    
    tag { params.prefix }

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}"

    input:
    file(tree_newick)


    output:
    path "*.nwk", emit: ch_reroot_iqtree

    script:
    """
    nw_reroot ${tree_newick} `head -1 ${params.reference_fasta} | tr -d \">\"` > reroot_phylogenetic_tree.nwk
    """
}

process MAKEALLELES {
    
    tag { params.prefix }

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}"

    input:
    file(msa_mafft)

    output:
    path "*.tsv", emit: ch_alleles

    script:
    """
    align2alleles.py --reference-name ${params.reference_name} ${msa_mafft} > alleles.tsv
    """
}

process VISUALIZE_PHYLOGENTICTREE {

    tag { params.prefix }

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}"

    input:
    file(newick_tree)
    file(lineage_report)
    file(alleles_info)

    output:
    path "*.pdf", emit: ch_visualization_phylogenetic_tree

    script:
    """
    phylogenetic_ggtree.r -n ${newick_tree} -l ${lineage_report} -a ${alleles_info}
    """
}

process VISUALIZE_SHIPTV_PHYLOGENTICTREE {

    tag { params.prefix }

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}"

    input:
    file(newick_tree)

    output:
    path "*.html", emit: ch_shiptv_html_visualization
    path "*.tsv", emit: ch_shiptv_metadata

    script:
    """
    shiptv -n ${newick_tree} -o shiptv_phylogenetic_tree.html -m shiptv_metadata.tsv
    """
}


