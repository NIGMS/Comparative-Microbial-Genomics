#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

def checkPathParamList = [ params.input, params.blast_db, params.busco_db, params.bakta_db ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

/*
========================================================================================
    IMPORT MODULES/SUBWORKFLOWS
========================================================================================
*/

include { INPUT_CHECK                 } from '../subworkflows/local/inputcheck/main'
include { FASTQC                      } from '../modules/nf-core/fastqc/main'
include { FASTP                       } from '../modules/nf-core/fastp/main'
include { QUAST                       } from '../modules/nf-core/quast/main'   
include { BWA_MEM                     } from '../modules/nf-core/bwa/mem/main'
include { BWA_INDEX                   } from '../modules/nf-core/bwa/index/main'
include { BUSCO                       } from '../modules/nf-core/busco/main'
include { SAMTOOLS_INDEX              } from '../modules/nf-core/samtools/index/main'
include { BAKTADBDOWNLOAD             } from '../modules/nf-core/bakta/baktadbdownload/main'
include { BAKTA                       } from '../modules/nf-core/bakta/bakta/main'
include { BLAST_BLASTN                } from '../modules/nf-core/blast/blastn/main'
include { SPADES                      } from '../modules/nf-core/spades/main'
include { BLOBTOOLS                   } from '../modules/local/blobtools/main'
include { paramsSummaryMap            } from 'plugin/nf-schema'

// Local modules
include { FINALIZE_RESULTS            } from '../modules/local/finalizeresults'

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

workflow WGSBAC {

    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        file(params.input)
    )
    ch_reads = INPUT_CHECK.out.reads
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    //
    // MODULE: Run FastQC
    //
    FASTQC(ch_reads)
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    //
    // MODULE: Run Fastp
    //
    FASTP(
        ch_reads,
        params.adapter_fasta,
        params.discard_trimmed_pass,
        params.save_trimmed_fail,
        params.save_merged
    )
    ch_versions = ch_versions.mix(FASTP.out.versions.first())

    //
    // MODULE: Run SPAdes
    //
    spades_input = FASTP.out.reads.map { meta, reads ->
        [meta, reads, [], []]  // Add empty lists for pacbio and nanopore
    }
    SPADES(
        spades_input,
        params.hmm,
        params.yml
    )
    ch_versions = ch_versions.mix(SPADES.out.versions.first())

    //
    // MODULE: Run BWA Index
    //
    BWA_INDEX(SPADES.out.contigs)
    ch_versions = ch_versions.mix(BWA_INDEX.out.versions.first())

    //
    // MODULE: Run BWA MEM
    //
    BWA_MEM(
        FASTP.out.reads,
        BWA_INDEX.out.index,
        SPADES.out.contigs,
        true
    )
    ch_versions = ch_versions.mix(BWA_MEM.out.versions.first())

    //
    // MODULE: Run Samtools Index
    //
    SAMTOOLS_INDEX(BWA_MEM.out.bam)
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    //
    // MODULE: Run BUSCO
    //
    busco_db_ch = Channel.fromPath(params.busco_db, checkIfExists: true)

    BUSCO(
        SPADES.out.contigs,
        params.mode,
        params.lineage,
        busco_db_ch,
        params.config_file,
        params.clean_intermediates)
    ch_versions = ch_versions.mix(BUSCO.out.versions.first())

    //
    // MODULE: Run QUAST
    //

    QUAST(
        SPADES.out.contigs
    )
    ch_versions = ch_versions.mix(QUAST.out.versions.first())

    //
    // MODULE: Run BLAST
    //

    blast_db_ch = Channel.of([[:], file(params.blast_db)])

    BLAST_BLASTN(
        SPADES.out.contigs,
        blast_db_ch
    )
    ch_versions = ch_versions.mix(BLAST_BLASTN.out.versions.first())

    //
    // MODULE: Run Blobtools
    //

    // Combine channels into a single channel for BLOBTOOLS
    
    ch_blobtools_input = SPADES.out.contigs.join(BWA_MEM.out.bam).join(SAMTOOLS_INDEX.out.bai).join(BLAST_BLASTN.out.txt)
    nodes_db_ch = Channel.of(file(params.nodes_db))

    BLOBTOOLS(
        ch_blobtools_input,
        nodes_db_ch
    )

    ch_versions = ch_versions.mix(BLOBTOOLS.out.versions.first())

    //
    // MODULE: Run Bakta
    //

    BAKTADBDOWNLOAD()

    BAKTA(
        SPADES.out.contigs,
        BAKTADBDOWNLOAD.out.db,
        params.bakta_proteins,
        params.bakta_prodigal_tf
    )
    ch_versions = ch_versions.mix(BAKTA.out.versions.first())

    //
    // MODULE: Finalize Results
    //

    // Collect outputs for FINALIZE_RESULTS
    ch_final_input = FASTQC.out.zip
        .join(SPADES.out.contigs)
        .join(QUAST.out.results)
        .join(BUSCO.out.short_summaries_txt)
        .join(BAKTA.out.faa)
        .join(BLOBTOOLS.out.plot)
        .join(BLOBTOOLS.out.cov_plot)


    FINALIZE_RESULTS(ch_final_input)
}

/*
========================================================================================
    COMPLETION EMAIL AND SUMMARY
========================================================================================
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log)
    }
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
========================================================================================
    THE END
========================================================================================
*/