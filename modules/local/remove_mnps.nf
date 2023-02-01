process REMOVE_MNPS{
    tag "${meta.id}"
    label 'process_medium'

    conda "bioconda::pysam=0.15.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'quay.io/biocontainers/pysam:v0.15.2ds-2-deb-py3_cv1':
        'quay.io/biocontainers/pysam:v0.15.2ds-2-deb-py3_cv1' }"

    input:
        tuple val(meta), path(vcf), path(idx)

    output:
    tuple val(meta), path("*.rmMnps.vcf.gz"), path("*.rmMnps.vcf.gz.tbi"),    emit: tuple_vcf_remove_mnps
    path "versions.yml",                                                        emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def avail_mem = 3
    if (!task.memory) {
        log.info '[GATK SelectVariants] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.toGiga()
    }
    """
    
    python3 $baseDir/bin/removeMnps.py -v ${vcf} -o ${meta.id}.rmMnps.vcf.gz
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pysam: \$(echo "0.15.2")
    END_VERSIONS
    """
}
