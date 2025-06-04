process BLOBTOOLS {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::blobtools=1.1.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/blobtools:1.1.1--py_1':
        'biocontainers/blobtools:1.1.1--py_1' }"

    input:
    tuple val(meta), path(genome), path(bam), path(bam_index), path(blast_results)
    path db

    output:
    tuple val(meta), path("*.blobDB.json")              , emit: blobdb
    tuple val(meta), path("*_table.tsv")                , emit: table
    tuple val(meta), path("*_blobplot.png")             , emit: plot
    tuple val(meta), path("*_blobplot_read_cov.png")    , emit: cov_plot
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Set matplotlib backend
    export MPLBACKEND=Agg

    # Decompress the genome file
    zcat ${genome} > ${prefix}.fa

    # Create database
    blobtools create \\
        -i ${prefix}.fa \\
        -b $bam \\
        -t $blast_results \\
        -o ${prefix} \\
        --db $db \\
        $args

    # Produce results table
    blobtools view \\
        -i ${prefix}.blobDB.json \\
        -r all \\
        -o ${prefix}

    # Generate figures
    blobtools plot \\
        -i ${prefix}.blobDB.json \\
        -r genus

    # Rename files to match expected output
    mv ${prefix}.${prefix}.blobDB.table.txt ${prefix}_table.tsv
    mv "${prefix}.blobDB.json.bestsum.genus.p8.span.100.blobplot.bam0.png" "${prefix}_blobplot.png"
    mv "${prefix}.blobDB.json.bestsum.genus.p8.span.100.blobplot.read_cov.bam0.png" "${prefix}_blobplot_read_cov.png"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blobtools: \$(blobtools --version 2>&1 | sed 's/blobtools v//')
        zcat: \$(zcat --version 2>&1 | sed 's/zcat (gzip) //')
    END_VERSIONS
    """
}