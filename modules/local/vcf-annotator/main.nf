process  VCF_ANNOTATOR{
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::vcf-annotator==0.7--hdfd78af_0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/vcf-annotator:0.7--hdfd78af_0':
        'biocontainers/vcf-annotator:0.7--hdfd78af_0' }"

    input:
    tuple val(meta), path(vcf)
    path  fasta
    path  genbank

    output:
    tuple val(meta), path("*.annotated.vcf")   , emit: vcf
    path "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def avail_mem = 3072
    if (!task.memory) {
        log.info '[GATK VariantsToTable] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }
    """
    python3 vcf-annotator.py \\
        ${vcf} \\
        ${fasta} \\
        ${genbank}
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vcf-annotator: \$(echo \$(vcf-annotator --version 2>&1) | sed 's/^.*(VCF_ANNOTATOR) v//; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.variants.table

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
    vcf-annotator: \$(echo \$(vcf-annotator --version 2>&1) | sed 's/^.*(VCF_ANNOTATOR) v//; s/ .*\$//')
    END_VERSIONS
    """
}

