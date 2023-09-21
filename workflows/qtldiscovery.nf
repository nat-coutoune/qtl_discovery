/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryLog; paramsSummaryMap } from 'plugin/nf-validation'

def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
def summary_params = paramsSummaryMap(workflow)

// Print parameter summary log to screen
log.info logo + paramsSummaryLog(workflow) + citation

WorkflowQtldiscovery.initialise(params, log)

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
include { INPUT_CHECK } from '../subworkflows/local/input_check'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC                          } from '../modules/nf-core/fastqc/main'
include { MULTIQC                         } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS     } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { TRIMMOMATIC                     } from '../modules/nf-core/trimmomatic/main'
include { PEAR                            } from '../modules/nf-core/pear/main'
include { BOWTIE2_BUILD                   } from '../modules/nf-core/bowtie2/build/main'
include { BOWTIE2_ALIGN                   } from '../modules/nf-core/bowtie2/align/main'
include { SAMTOOLS_SORT                   } from '../modules/nf-core/samtools/sort/main'
include { PICARD_MARKDUPLICATES           } from '../modules/nf-core/picard/markduplicates/main'                           
include { SAMTOOLS_FAIDX                  } from '../modules/nf-core/samtools/faidx/main'         
include { PICARD_CREATESEQUENCEDICTIONARY } from '../modules/nf-core/picard/createsequencedictionary/main'
include { SAMTOOLS_INDEX                  } from '../modules/nf-core/samtools/index/main'                                        
include { GATK4_HAPLOTYPECALLER           } from '../modules/nf-core/gatk4/haplotypecaller/main'      
include { GATK4_VARIANTFILTRATION         } from '../modules/nf-core/gatk4/variantfiltration/main'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow QTLDISCOVERY {

    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        file(params.input)
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)
    // TODO: OPTIONAL, you can use nf-validation plugin to create an input channel from the samplesheet with Channel.fromSamplesheet("input")
    // See the documentation https://nextflow-io.github.io/nf-validation/samplesheets/fromSamplesheet/
    // ! There is currently no tooling to help you write a sample sheet schema

    //
    // MODULE: Run FastQC
    //
    FASTQC (
        INPUT_CHECK.out.reads
    )
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: TRIMMOMATIC
    //
    TRIMMOMATIC(
      INPUT_CHECK.out.reads
    )
    ch_versions = ch_versions.mix(TRIMMOMATIC.out.versions.first())
    
    //
    // MODULE: PEAR
    //
    TRIMMOMATIC.out.trimmed_reads
      .branch { meta, reads ->
        pe: meta.single_end == false
        se: true
      }
      .set { all_reads }
  
    PEAR(
      all_reads.pe
    )
    
    PEAR.out.assembled
      .mix(all_reads.se)
      .map { meta, read -> [[id:meta.id, single_end:true], read] }
      .set { consensus_reads }
    
    //return

    ch_versions = ch_versions.mix(PEAR.out.versions.first())

    //
    // MODULE: BOWTIE2_BUILD
    //
    Channel
      .fromPath(params.fasta, checkIfExists: true)
      .map { it -> [[:], it] } //dicionario para dar formato para do bowtie2-build
      .set { fasta }

    BOWTIE2_BUILD(
      fasta
    )
    ch_versions = ch_versions.mix(BOWTIE2_BUILD.out.versions.first())

    //
    // MODULE: BOWTIE2 ALIGNER
    //
    index = BOWTIE2_BUILD.out.index.first()

    BOWTIE2_ALIGN(
      consensus_reads,
      index,
      true,
      true
    )
    ch_versions = ch_versions.mix(BOWTIE2_ALIGN.out.versions.first())

    //
    // MODULE: SAMTOOLS FASTA INDEX
    //
    SAMTOOLS_FAIDX(
      fasta.first(),
      [[],[]]
    )
    ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions.first())

    //
    // MODULE: PICARD MARK DUPLICATES
    //
    bam   = BOWTIE2_ALIGN.out.aligned
    fai   = SAMTOOLS_FAIDX.out.fai

    PICARD_MARKDUPLICATES(
      bam,
      fasta.first(),
      fai.first()
     )
    ch_versions = ch_versions.mix(PICARD_MARKDUPLICATES.out.versions.first())

    //
    // MODULE: PICARD_CREATESEQUENCEDICTIONARY
    //
    Channel
      .fromPath(params.fasta, checkIfExists: true)
      .map { it -> [[id: 'genome'], it] } 
      .set { fasta }

    PICARD_CREATESEQUENCEDICTIONARY(
      fasta
    )
    ch_versions = ch_versions.mix(PICARD_CREATESEQUENCEDICTIONARY.out.versions.first())

    //
    // MODULE: SAMTOOLS_INDEX
    //
    bam = BOWTIE2_ALIGN.out.aligned

    SAMTOOLS_INDEX(
      bam
    )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())


    //
    // MODULE: GATK4_HAPLOTYPECALLER
    //
    BOWTIE2_ALIGN.out.aligned
      .mix(SAMTOOLS_INDEX.out.bai)
      .groupTuple()
      .map { meta, files -> [meta, files[0], files[1], [], []] }
      .set { meta_inp_inpindex }

    fai = SAMTOOLS_FAIDX.out.fai
    dict = PICARD_CREATESEQUENCEDICTIONARY.out.reference_dict

    GATK4_HAPLOTYPECALLER(
      meta_inp_inpindex,
      fasta.map { meta, fasta -> fasta }.first(),
      fai.map { meta, fai -> fai }.first(),
      dict.map { meta, dict -> dict }.first(),
      [],
      []
    )
    ch_versions = ch_versions.mix(GATK4_HAPLOTYPECALLER.out.versions.first())

    //
    // MODULE: GATK4_VARIANTFILTRATION
    //
    GATK4_HAPLOTYPECALLER.out.vcf
       .mix(GATK4_HAPLOTYPECALLER.out.tbi)
       .groupTuple()
       .map { meta, files -> [meta, files[0], files[1]] }
       .set { meta_vcf_filtration }

     Channel
      .fromPath(params.fasta, checkIfExists: true)
      .map { it -> [[:], it] } 
      .set { fasta }

     SAMTOOLS_FAIDX.out.fai
       .map { it -> [[:], it] } 
       .set { fai }

    PICARD_CREATESEQUENCEDICTIONARY.out.reference_dict
      .map { it -> [[:], it] } 
      .set { dict }
 
    GATK4_VARIANTFILTRATION ( 
       meta_vcf_filtration,
       fasta,
       fai,
       dict
     )
     ch_versions = ch_versions.mix(GATK4_VARIANTFILTRATION.out.versions.first())
/*     return

    //
    // MODULE: GATK4_COMBINEGVCFS 
    //
    GATK4_VARIANTFILTRATION.out.vcf
         .mix(GATK4_VARIANTFILTRATION.out.tbi)
         .groupTuple()
         .map { meta, files -> [meta, files[0], files[1]] }
         .set { meta_vcf_combine }

     GATK4_COMBINEGVCFS (
       meta_vcf_combine,
       fasta.map { meta, fasta -> fasta }.first(),
       fai.map { meta, fai -> fai }.first(),
       dict.map { meta, dict -> dict }.first()
     )
     
    //
    // MODULE: GATK4_VARIANTSTOTABLE
    //
    GATK4_HAPLOTYPECALLER.out.vcf
      .mix(GATK4_HAPLOTYPECALLER.out.tbi)
      .groupTuple()
      .map { meta, files -> [meta, files[0], files[1]] }
      .set { meta_vcf_totable }


    GATK4_VARIANTSTOTABLE(
      meta_vcf_totable,
      fasta.map { meta, fasta -> fasta }.first(),
      fai.map { meta, fai -> fai }.first(),
      dict.map { meta, dict -> dict }.first()
    )
    ch_versions = ch_versions.mix(GATK4_VARIANTSTOTABLE.out.versions.first())

*/
    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowQtldiscovery.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowQtldiscovery.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description, params)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )
    multiqc_report = MULTIQC.out.report.toList()
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
