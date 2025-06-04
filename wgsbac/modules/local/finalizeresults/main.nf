process FINALIZE_RESULTS {
    tag "$meta.id"
    publishDir "${params.outdir}/final_reports", mode: 'copy'

    input:
    tuple val(meta), path(fastqc_results), path(contigs), path(quast_results), path(busco_results), path(bakta_results), path(blob_plot), path(blob_cov_plot)

    output:
    tuple val(meta), path("${meta.id}_final_report.txt"), emit: report
    path "versions.yml"                                 , emit: versions
    path "final_results"                                , emit: all_results

    script:
    """
    mkdir -p final_results/fastqc_results
    mkdir -p final_results/contigs
    mkdir -p final_results/busco_results
    mkdir -p final_results/quast_results
    mkdir -p final_results/bakta_results/proteomes
    mkdir -p final_results/blobtools_plots

    echo "Sample ID: ${meta.id}" > ${meta.id}_final_report.txt
    echo "FASTQC Results: final_results/fastqc_results" >> ${meta.id}_final_report.txt
    echo "Contigs: final_results/contigs" >> ${meta.id}_final_report.txt
    echo "QUAST Results: final_results/quast_results" >> ${meta.id}_final_report.txt
    echo "BUSCO Results: final_results/busco_results" >> ${meta.id}_final_report.txt
    echo "BAKTA Results: final_results/bakta_results/proteomes" >> ${meta.id}_final_report.txt
    echo "Blobtools Plot: final_results/blobtools_plots" >> ${meta.id}_final_report.txt
    echo "Blobtools COV Plot: final_results/blobtools_plots" >> ${meta.id}_final_report.txt

    cp -R $fastqc_results final_results/fastqc_results/
    cp $contigs final_results/contigs/
    cp -R $quast_results final_results/quast_results/
    cp -R $busco_results final_results/busco_results/
    cp $bakta_results final_results/bakta_results/proteomes/
    cp $blob_plot final_results/blobtools_plots/
    cp $blob_cov_plot final_results/blobtools_plots/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        finalize_results: \$(echo "1.0" | sed 's/^/v/')
    END_VERSIONS
    """
}