process MERGE_CONSENSUS_SEQUENCES {

    tag { params.prefix }

    output:
    path "*.fasta", emit: ch_mergeconsensus
    
    script:
    """
    process_consensus.py --identifier ${params.consensus_identifier} --consensus_path ${params.consensus_path} --reference_fasta ${params.reference_fasta} --merge_consensus_sequence merge_consensus_sequence.fasta
    """
}

process MSA_MAFFT {
    
    tag { params.prefix }

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}"

    input:
    file(cat_consensus_seqeuences)

    output:
    path "msa.aln", emit: ch_msa_mafft

    script:
    """
    mafft \\
        --thread ${task.cpus} \\
        --6merpair \\
        --addfragments ${cat_consensus_seqeuences} \\
        ${params.reference_fasta} > msa.aln
    """
}
