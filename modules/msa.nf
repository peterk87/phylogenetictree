process catConsensusSequences {

    tag { params.prefix }

    input:
    path (samples_consensus)

    output:
    path "*.fasta", emit: ch_catconsensus
    
    script:
    """
    for files in ${samples_consensus}
    do 
        cat \$files >> mergeConsensusSequences.fasta
    done
    """
}

process msaMAFFT {
    
    tag { params.prefix }

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}"

    input:
    file(cat_consensus_seqeuences)

    output:
    path "msa.aln", emit: ch_msa_mafft

    script:
    """
    mafft --ep --localpair  --maxiterate 16 --reorder ${cat_consensus_seqeuences} > "msa.aln"
    """
}




