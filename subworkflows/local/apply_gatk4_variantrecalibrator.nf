//
// Check input samplesheet and get read channels
//

include { GATK4_VARIANTRECALIBRATOR } from '../../modules/nf-core/gatk4/variantrecalibrator/main'

workflow APPLY_GATK4_VARIANTRECALIBRATOR{
    take:
        vcf_ref_fai_seq

    main:
        Channel
            .fromPath( params.resource_vcfs_file )
            .splitText()
            .map{ it -> file( it.trim() ) }
            .collect()
            .set{ resource_vcfs }
        Channel
            .fromPath( params.resource_vcf_idx_file )
            .splitText()
            .map{ it -> file( it.trim() ) }
            .collect()
            .set { resource_vcf_idx }
        Channel
            .fromPath( params.resource_labels_file )
            .splitText()
            .map{ it -> it.trim() }
            .collect()
            .set { resource_labels }

       tuple1 = vcf_ref_fai_seq.combine( resource_vcfs )
       tuple2 = tuple1.combine( resource_vcf_idx )
       tuple3 = tuple2.combine( resource_labels )

       tuple3.view()

        vcf_input     = tuple3.map{ meta, vcf, idx, ref, fai, seq, resVcfs, resIdx, resLabel -> tuple(meta, vcf, idx) }
        //resVcf_input  = tuple3.map{ vcf, ref, fai, seq, resVcfs, resIdx, resLabel -> resVcfs }
        //resTbi_input  = tuple3.map{ vcf, ref, fai, seq, resVcfs, resIdx, resLabel -> resIdx }
        //resLab_input  = tuple3.map{ vcf, ref, fai, seq, resVcfs, resIdx, resLabel -> resLabel }
        ref_input     = tuple3.map{ meta, vcf, idx, ref, fai, seq, resVcfs, resIdx, resLabel -> ref }
        fai_input     = tuple3.map{ meta, vcf, idx, ref, fai, seq, resVcfs, resIdx, resLabel -> fai }
        seqDict_input = tuple3.map{ meta, vcf, idx, ref, fai, seq, resVcfs, resIdx, resLabel -> seq }

        GATK4_VARIANTRECALIBRATOR(
         vcf_input,
         resource_vcfs,
         resource_vcf_idx,
         resource_labels,
         ref_input,
         fai_input,
         seqDict_input
        )
}

