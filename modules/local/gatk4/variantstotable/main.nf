process GATK4_VARIANTSTOTABLE {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::gatk4==4.4.0.0--py36hdfd78af_0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.4.0.0--py36hdfd78af_0' :
        'biocontainers/gatk4:4.4.0.0--py36hdfd78af_0' }"


    input:
    tuple val(meta), path(vcf), path (tbi)
    path  fasta
    path  fai
    path  dict


    output:
    tuple val(meta), path("*.variants.table")   , emit: tsv
    path "versions.yml"                         , emit: versions

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
    gatk VariantsToTable \\
        -V ${vcf} \\
        -F CHROM -F POS -F TYPE -F FILTER -GF DP -GF AD \\
        -R ${fasta} \\
        -O ${prefix}.variants.table
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
/*
    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def bamout_command = args.contains("--bam-writer-type") ? "--bam-output ${prefix.replaceAll('.g\\s*$', '')}vcf.gz" : ""

    def stub_realigned_bam = bamout_command ? "touch ${prefix.replaceAll('.g\\s*$', '')}.vcf.gz" : ""
    """
    touch ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
*/
}

