#!/usr/bin/env nextflow

/*
========================================================================================
    PLASMID ANALYSIS PIPELINE (PlasmidScope)
========================================================================================
    A comprehensive Nextflow pipeline for bacterial plasmid characterization,
    mobile genetic element analysis, and genomic environment visualization.

    GitHub  : https://github.com/<username>/plasmid-analysis-pipeline
    Docs    : docs/
    License : MIT
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

// ========================================================================================
//    PRINT PIPELINE HEADER
// ========================================================================================

log.info """\
╔═══════════════════════════════════════════════════════════════════╗
║                                                                   ║
║   ██████╗ ██╗      █████╗ ███████╗███╗   ███╗██╗██████╗          ║
║   ██╔══██╗██║     ██╔══██╗██╔════╝████╗ ████║██║██╔══██╗         ║
║   ██████╔╝██║     ███████║███████╗██╔████╔██║██║██║  ██║         ║
║   ██╔═══╝ ██║     ██╔══██║╚════██║██║╚██╔╝██║██║██║  ██║         ║
║   ██║     ███████╗██║  ██║███████║██║ ╚═╝ ██║██║██████╔╝         ║
║   ╚═╝     ╚══════╝╚═╝  ╚═╝╚══════╝╚═╝     ╚═╝╚═╝╚═════╝          ║
║                                                                   ║
║   ███████╗ ██████╗ ██████╗ ██████╗ ███████╗                      ║
║   ██╔════╝██╔════╝██╔═══██╗██╔══██╗██╔════╝                      ║
║   ███████╗██║     ██║   ██║██████╔╝█████╗                        ║
║   ╚════██║██║     ██║   ██║██╔═══╝ ██╔══╝                        ║
║   ███████║╚██████╗╚██████╔╝██║     ███████╗                      ║
║   ╚══════╝ ╚═════╝ ╚═════╝ ╚═╝     ╚══════╝                      ║
║                                                                   ║
║   PlasmidScope v${workflow.manifest.version}                      ║
║   Comprehensive Plasmid Genomic Analysis Pipeline                 ║
║                                                                   ║
╚═══════════════════════════════════════════════════════════════════╝

    Pipeline Parameters:
    ====================
    Input samplesheet : ${params.input}
    Output directory  : ${params.outdir}
    Assembly mode     : ${params.assembly_mode}
    Reference genome  : ${params.reference ?: 'Not provided'}
    Skip assembly     : ${params.skip_assembly}
    Skip MGE analysis : ${params.skip_mge}
    Skip visualization: ${params.skip_viz}
    Max CPUs          : ${params.max_cpus}
    Max Memory        : ${params.max_memory}
    """.stripIndent()

// ========================================================================================
//    VALIDATE INPUTS
// ========================================================================================

if (params.input) {
    ch_input = Channel.fromPath(params.input, checkIfExists: true)
} else {
    exit 1, "ERROR: Please provide an input samplesheet with --input"
}

// ========================================================================================
//    IMPORT MODULES
// ========================================================================================

include { INPUT_CHECK           } from './modules/input_check'
include { FASTQC                } from './modules/fastqc'
include { NANOPLOT              } from './modules/nanoplot'
include { FASTP                 } from './modules/fastp'
include { FILTLONG              } from './modules/filtlong'
include { UNICYCLER             } from './modules/unicycler'
include { FLYE                  } from './modules/flye'
include { HYBRACTER             } from './modules/hybracter'
include { QUAST                 } from './modules/quast'
include { BAKTA                 } from './modules/bakta'
include { ABRICATE              } from './modules/abricate'
include { AMRFINDERPLUS         } from './modules/amrfinderplus'
include { PLASMIDFINDER         } from './modules/plasmidfinder'
include { MOB_SUITE             } from './modules/mob_suite'
include { ISFINDER_BLAST        } from './modules/isfinder_blast'
include { INTEGRON_FINDER       } from './modules/integron_finder'
include { MOBILE_ELEMENT_FINDER } from './modules/mobile_element_finder'
include { BRIG                  } from './modules/brig'
include { CLINKER               } from './modules/clinker'
include { MASH_DIST             } from './modules/mash_dist'
include { PLASMID_REPORT        } from './modules/plasmid_report'
include { MULTIQC               } from './modules/multiqc'

// ========================================================================================
//    NAMED WORKFLOWS
// ========================================================================================

workflow PLASMIDSCOPE {

    take:
    samplesheet

    main:

    // -------------------------------------------------------
    // PHASE 1: Input Validation and QC
    // -------------------------------------------------------
    INPUT_CHECK ( samplesheet )
    ch_reads = INPUT_CHECK.out.reads

    // Split into short and long reads
    ch_short_reads = ch_reads.map { meta, sr1, sr2, lr ->
        [ meta, sr1, sr2 ]
    }.filter { meta, sr1, sr2 -> sr1 != null }

    ch_long_reads = ch_reads.map { meta, sr1, sr2, lr ->
        [ meta, lr ]
    }.filter { meta, lr -> lr != null }

    // Short-read QC
    FASTQC ( ch_short_reads )
    FASTP ( ch_short_reads )

    // Long-read QC
    NANOPLOT ( ch_long_reads )
    FILTLONG ( ch_long_reads )

    // -------------------------------------------------------
    // PHASE 2: Genome Assembly
    // -------------------------------------------------------
    if (!params.skip_assembly) {

        if (params.assembly_mode == 'hybrid') {
            // Combine trimmed short reads with filtered long reads
            ch_assembly_input = FASTP.out.reads
                .join( FILTLONG.out.reads )

            UNICYCLER ( ch_assembly_input )
            ch_assembly = UNICYCLER.out.assembly

        } else if (params.assembly_mode == 'long') {
            FLYE ( FILTLONG.out.reads )
            ch_assembly = FLYE.out.assembly

        } else if (params.assembly_mode == 'auto') {
            HYBRACTER ( FASTP.out.reads.join(FILTLONG.out.reads) )
            ch_assembly = HYBRACTER.out.assembly
        }

        QUAST ( ch_assembly )

    } else {
        // Use pre-assembled genomes from samplesheet
        ch_assembly = INPUT_CHECK.out.assemblies
    }

    // -------------------------------------------------------
    // PHASE 3: Annotation
    // -------------------------------------------------------
    BAKTA ( ch_assembly )

    // -------------------------------------------------------
    // PHASE 4: AMR and Virulence Screening
    // -------------------------------------------------------
    ABRICATE ( ch_assembly )
    AMRFINDERPLUS ( ch_assembly, BAKTA.out.proteins, BAKTA.out.gff )

    // -------------------------------------------------------
    // PHASE 5: Plasmid Typing
    // -------------------------------------------------------
    PLASMIDFINDER ( ch_assembly )
    MOB_SUITE ( ch_assembly )

    // -------------------------------------------------------
    // PHASE 6: Mobile Genetic Element Analysis
    // -------------------------------------------------------
    if (!params.skip_mge) {
        ISFINDER_BLAST ( ch_assembly )
        INTEGRON_FINDER ( ch_assembly )
        MOBILE_ELEMENT_FINDER ( ch_assembly )
    }

    // -------------------------------------------------------
    // PHASE 7: Comparative Analysis and Visualization
    // -------------------------------------------------------
    if (!params.skip_viz) {
        // Collect all plasmid sequences for comparison
        ch_plasmids = MOB_SUITE.out.plasmids.collect()

        MASH_DIST ( ch_plasmids )

        if (params.reference) {
            ch_reference = Channel.fromPath(params.reference)
            BRIG ( ch_reference, ch_plasmids )
        }

        CLINKER ( BAKTA.out.gbk.collect() )
    }

    // -------------------------------------------------------
    // PHASE 8: Reporting
    // -------------------------------------------------------
    PLASMID_REPORT (
        BAKTA.out.summary.collect(),
        ABRICATE.out.results.collect(),
        AMRFINDERPLUS.out.results.collect(),
        PLASMIDFINDER.out.results.collect(),
        MOB_SUITE.out.results.collect(),
        params.skip_mge ? Channel.empty() : ISFINDER_BLAST.out.results.collect(),
        params.skip_mge ? Channel.empty() : INTEGRON_FINDER.out.results.collect()
    )

    MULTIQC (
        FASTQC.out.zip.collect().ifEmpty([]),
        FASTP.out.json.collect().ifEmpty([]),
        QUAST.out.results.collect().ifEmpty([])
    )

    emit:
    assembly    = ch_assembly
    report      = PLASMID_REPORT.out.report
    multiqc     = MULTIQC.out.report
}

// ========================================================================================
//    RUN MAIN WORKFLOW
// ========================================================================================

workflow {
    PLASMIDSCOPE ( ch_input )
}

// ========================================================================================
//    ON COMPLETE
// ========================================================================================

workflow.onComplete {
    log.info """\
    ╔═══════════════════════════════════════════════════════╗
    ║   PlasmidScope Pipeline Complete!                     ║
    ╠═══════════════════════════════════════════════════════╣
    ║   Status    : ${workflow.success ? 'SUCCESS' : 'FAILED'}
    ║   Duration  : ${workflow.duration}
    ║   Output    : ${params.outdir}
    ╚═══════════════════════════════════════════════════════╝
    """.stripIndent()
}
