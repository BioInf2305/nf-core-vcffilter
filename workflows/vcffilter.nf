/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
//WorkflowFiltervcf.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
//def checkPathParamList = [ params.input, params.multiqc_config, params.fasta ]
//for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
//if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { FASTA_INDICES } from '../subworkflows/local/fasta_indices'
include { APPLY_GATK4_VARIANTRECALIBRATOR } from '../subworkflows/local/apply_gatk4_variantrecalibrator.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//

include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { GATK4_SELECTVARIANT } from '../modules/local/gatk4/selectvariants/main'
include { REMOVE_MNPS } from '../modules/local/remove_mnps'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow VCFFILTER {

    ch_versions = Channel.empty()

    Channel
        .fromFilePairs(params.input)
        .view()
        .set{ vcfPair }
    
    // convert it into meta format: [meta, vcf, idx]

    tuple_meta_vcf = vcfPair.map{ sampleN, vcfPa -> tuple([id:sampleN, single_end:false], vcfPa[0], vcfPa[1]) }

    //
    // module: gatk4_selectvariants 
    //

    if ( !params.skip_gatk4_selectvariants ){
        GATK4_SELECTVARIANT(tuple_meta_vcf)
        out1 = GATK4_SELECTVARIANT.out.tuple_vcf_gatk4_selectvariant
        ch_versions = ch_versions.mix(GATK4_SELECTVARIANT.out.versions)
    }
    else{
         out1= tuple_meta_vcf
        }

    //
    // module: remove_mnps
    //

    if ( !params.skip_remove_mnps ){
        REMOVE_MNPS(out1)
        out2 = REMOVE_MNPS.out.tuple_vcf_remove_mnps
        ch_versions = ch_versions.mix(REMOVE_MNPS.out.versions)
    }
    else{
        out2 = out1
        }

    //
    // SUBWORKFLOW: Deals with bwa index and samtools faidx
    //
    
    fastaF = Channel.fromPath(params.fasta)

    FASTA_INDICES(
        fastaF
     )
    
    picard_ver = params.skip_picard_seqdict ?'':FASTA_INDICES.out.picard_seqdict_version

    faidx_ver = params.skip_samtools_faidx?'':FASTA_INDICES.out.samtools_faidx_version

    if ( picard_ver != ''){
        ch_versions = ch_versions.mix( picard_ver )
        }
    if ( faidx_ver != ''){
        ch_version = ch_versions.mix(faidx_ver)
        }
    
    // prepare input channel for GATK4_VARIANTRECALIBRATOR

    vcf_fasta = out2.combine( fastaF )
    vcf_fasta_fai = vcf_fasta.combine( FASTA_INDICES.out.fa_idx )
    vcf_fasta_fai_seqdict = vcf_fasta_fai.combine( FASTA_INDICES.out.picard_idx )


    //
    // SUBWORKFLOW: APPLY_GATK4_VARIANTRECALIBRATOR
    //

    APPLY_GATK4_VARIANTRECALIBRATOR(
        vcf_fasta_fai_seqdict
    )

    //
    // MODULE: MultiQC
    //
    //workflow_summary    = WorkflowFiltervcf.paramsSummaryMultiqc(workflow, summary_params)
    //ch_workflow_summary = Channel.value(workflow_summary)

    //methods_description    = WorkflowFiltervcf.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description)
    //ch_methods_description = Channel.value(methods_description)
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
