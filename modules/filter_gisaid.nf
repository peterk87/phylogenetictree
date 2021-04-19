
process FILTER_GISIAD_SEQUENCES {
    
    tag { params.prefix }

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}"

    input:
    file(ch_gisaid_sequence)
    file(ch_gisaid_metadata)

    output:
    path "*.fasta", emit: ch_filter_gisiad_sequences
    path "*.tsv", emit: ch_filter_gisiad_metadata


    script:
    filtered_fasta_output ="filtered_gisiad_sequences.fasta"
    filtered_metadata_output ="filtered_metadata.tsv"
    """
    filter_gisaid_sequences.py -i $ch_gisaid_sequence -m $ch_gisaid_metadata -of $filtered_fasta_output -om $filtered_metadata_output
    """
}
