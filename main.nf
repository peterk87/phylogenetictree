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

log.info nfcoreHeader()

def helpMessage() {
    // TODO nf-core: Add to this help message with new command line parameters
    
    log.info"""
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run nf-core/phylogenetictree -profile docker
    """.stripIndent()
}

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

custom_runName = params.name
if (!(workflow.runName ==~ /[a-z]+_[a-z]+/)) {
    custom_runName = workflow.runName
}

def summary = [:]
if (workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Run Name']         = custom_runName ?: workflow.runName
// TODO nf-core: Report custom parameters here
summary['Input']                            = params.input
summary['Reference Name']                   = params.reference_name
summary['Fasta Ref']                        = params.reference_fasta
summary['Max Resources']                    = "$params.max_memory memory, $params.max_cpus cpus, $params.max_time time per job"
if (workflow.containerEngine) summary['Container'] = "$workflow.containerEngine - $workflow.container"
summary['Output dir']       = params.outdir
summary['Launch dir']       = workflow.launchDir
summary['Working dir']      = workflow.workDir
summary['Script dir']       = workflow.projectDir
summary['User']             = workflow.userName
if (workflow.profile.contains('awsbatch')) {
    summary['AWS Region']   = params.awsregion
    summary['AWS Queue']    = params.awsqueue
    summary['AWS CLI']      = params.awscli
}
summary['Config Profile'] = workflow.profile
if (params.config_profile_description) summary['Config Profile Description'] = params.config_profile_description
if (params.config_profile_contact)     summary['Config Profile Contact']     = params.config_profile_contact
if (params.config_profile_url)         summary['Config Profile URL']         = params.config_profile_url
summary['Config Files'] = workflow.configFiles.join(', ')
if (params.email || params.email_on_fail) {
    summary['E-mail Address']    = params.email
    summary['E-mail on failure'] = params.email_on_fail
    summary['MultiQC maxsize']   = params.max_multiqc_email_size
}
log.info summary.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")
log.info "-\033[2m--------------------------------------------------\033[0m-"

// Check the hostnames against configured profiles
// checkHostname()

Channel.from(summary.collect{ [it.key, it.value] })
    .map { k,v -> "<dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }
    .reduce { a, b -> return [a, b].join("\n            ") }
    .map { x -> """
    id: 'nf-core-phylogenetictree-summary'
    description: " - this information is collected when the pipeline is started."
    section_name: 'nf-core/phylogenetictree Workflow Summary'
    section_href: 'https://github.com/nf-core/phylogenetictree'
    plot_type: 'html'
    data: |
        <dl class=\"dl-horizontal\">
            $x
        </dl>
    """.stripIndent() }
    .set { ch_workflow_summary }

include {FILTER_GISIAD_SEQUENCES} from './modules/filter_gisaid'
include {PHYLOGENETIC}  from './phylogenetic.nf'


workflow {
  
    ch_consensus_seqs = Channel
        .fromPath(params.input)
        .splitFasta( record: [id: true, sequence: true])
        .collectFile( name: 'consensus_seqs.fa' ){
        ">${it.id}\n${it.sequence}"
    }

    ch_reference_fasta = Channel.fromPath(params.reference_fasta)

    if (params.filter_gisaid){

        ch_gisaid_sequence = Channel.fromPath(params.gisiad_sequences)
    
        ch_gisaid_metadata = Channel.fromPath(params.gisiad_metadata)

        FILTER_GISIAD_SEQUENCES(ch_gisaid_sequence, ch_gisaid_metadata)

        PHYLOGENETIC(FILTER_GISIAD_SEQUENCES.out.ch_filter_gisiad_sequences, ch_reference_fasta)
    } 
    else{

        PHYLOGENETIC(ch_consensus_seqs, ch_reference_fasta)
    }

    
}


def nfcoreHeader() {
    // Log colors ANSI codes
    c_black = params.monochrome_logs ? '' : "\033[0;30m";
    c_blue = params.monochrome_logs ? '' : "\033[0;34m";
    c_cyan = params.monochrome_logs ? '' : "\033[0;36m";
    c_dim = params.monochrome_logs ? '' : "\033[2m";
    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_reset = params.monochrome_logs ? '' : "\033[0m";
    c_white = params.monochrome_logs ? '' : "\033[0;37m";
    c_yellow = params.monochrome_logs ? '' : "\033[0;33m";
    return """    -${c_dim}--------------------------------------------------${c_reset}-
                                            ${c_green},--.${c_black}/${c_green},-.${c_reset}
    ${c_blue}        ___     __   __   __   ___     ${c_green}/,-._.--~\'${c_reset}
    ${c_blue}  |\\ | |__  __ /  ` /  \\ |__) |__         ${c_yellow}}  {${c_reset}
    ${c_blue}  | \\| |       \\__, \\__/ |  \\ |___     ${c_green}\\`-._,-`-,${c_reset}
                                            ${c_green}`._,._,\'${c_reset}
    ${c_purple}  nf-core/phylogenetictree v${workflow.manifest.version}${c_reset}
    -${c_dim}--------------------------------------------------${c_reset}-
    """.stripIndent()
}