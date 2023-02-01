process GATK4_SELECTVARIANT {
    tag "${meta.id}"
    label 'process_medium'

    conda "bioconda::gatk4=4.3.0.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.3.0.0--py36hdfd78af_0':
        'quay.io/biocontainers/gatk4:4.3.0.0--py36hdfd78af_0' }"

    input:
        tuple val(meta), path(vcf), path(idx)

    output:
    tuple val(meta), path("*.selectvariants.vcf.gz"), path("*.selectvariants.vcf.gz.tbi"), emit: tuple_vcf_gatk4_selectvariant
    path "versions.yml",                                                                   emit: versions

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
    gatk --java-options "-Xmx${avail_mem}G" SelectVariants \\
        --variant $vcf \\
        --output ${meta.id}.selectvariants.vcf.gz \\
        --tmp-dir . \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    """
    touch ${prefix}.selectvariants.vcf.gz
    touch ${prefix}.selectvariants.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
