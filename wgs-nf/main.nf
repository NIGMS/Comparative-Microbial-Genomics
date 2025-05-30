#!/usr/bin/env nextflow

params.reads = './raw-reads/*_{1,2}.fastq.gz'
params.output_dir = './output-dir'
params.threads = 32
params.blast_db = "${launchDir}/databases/blast_db/blast_db"
// params.blast_db = "${launchDir}/databases/blast_db/ref_prok_rep_genomes"
params.busco_db = "${launchDir}/databases/busco_downloads/"
params.bakta_db = "${launchDir}/databases/bakta_db/"
params.blob_taxdb = "${launchDir}/databases/blast_db/nodesDB.txt"
params.blob_lookup = "${launchDir}/databases/blast_db/accessions_to_taxids.txt"

// Step 0: Count number of reads and determine read length
process ASSESS_READS {
    input:
    tuple val(sampleid), path(fq1), path(fq2)

    output:
    tuple val(sampleid), env('read_length'), env('num_reads'), env('throughput')


    script:
    """
    # Calculate the read length (example using first read in fq1)
    read_length=\$(zcat ${fq1} | head -n 2 | tail -n 1 | wc -c)

    # Calculate the number of reads (divide total lines by 4)
    num_reads=\$(zcat ${fq1} | wc -l | awk '{print \$1 / 4}')

    # Calculate the throughput for paired-end reads
    throughput=\$((read_length * num_reads * 2))

    # Output the results in a Nextflow-compatible format
    echo "${sampleid},\$read_length,\$num_reads,\$throughput"
    """
}



// Step 1: Quality control with FastQC
process RUN_FASTQC {
//     container 'quay.io/biocontainers/fastqc:0.11.9--hdfd78af_1'

    input:
    tuple val(sampleid), path(fq1), path(fq2)

    output:
    tuple val(sampleid), path( "output-fastqc/*.html" ), path( "output-fastqc/*.zip" )


    publishDir "${params.output_dir}/output-fastqc", pattern: "output-fastqc/*.html", mode: 'copy'

    script:
    """
    mkdir -p "output-fastqc/"
    
    fastqc ${fq1} ${fq2} \
	--outdir output-fastqc/ \
	--threads ${params.threads}
    """
}

// Step 2: Trim paired-end reads with FASTP
process RUN_FASTP {
    input:
    tuple val(sampleid), path(fq1), path(fq2)
    
    output:
    tuple val(sampleid), path("${sampleid}.trimmed_R1.fastq.gz"), path("${sampleid}.trimmed_R2.fastq.gz")
    
    script:
    """
    fastp -i $fq1 -I $fq2 \
        -o "$sampleid".trimmed_R1.fastq.gz -O "$sampleid".trimmed_R2.fastq.gz \
        --thread ${params.threads}
    """
}

// Step 3: Assemble genome with spades
process RUN_SPADES {
    errorStrategy 'ignore'

    input:
    tuple val(sampleid), path(fq1), path(fq2)

    output:
    tuple val(sampleid), path("${sampleid}.fasta")

    publishDir "${params.output_dir}/output-genome/", mode: 'copy', saveAs: { "${sampleid}.fasta" }

    script:
    """
    spades.py -1 $fq1 -2 $fq2 \
	-o output-spades-${sampleid} \
	--threads ${params.threads}

	mv output-spades-${sampleid}/contigs.fasta ${sampleid}.fasta


    """
}


process RUN_BWA {
    input:
    tuple val(sampleid), path(fq1), path(fq2), path(genome)

    output:
    tuple val(sampleid), path("${sampleid}.mapped.bam"), path("${sampleid}.mapped.bam.bai")
    
    script:
    """
    bwa index $genome

    bwa mem $genome $fq1 $fq2 \
	-t ${params.threads} \
        | samtools view -bS -@ ${params.threads} - \
	| samtools sort -@ ${params.threads} -o ${sampleid}.mapped.bam

    # index output file
    samtools index ${sampleid}.mapped.bam
    """
}

process RUN_BUSCO {
    input:
    tuple val(sampleid), path(genome)

    output:
    tuple val(sampleid), path("output-busco-${sampleid}/short_summary.specific.bacteria_odb10.output-busco-${sampleid}.txt")

    publishDir "${params.output_dir}/output-busco/", mode: 'copy', saveAs: { "${sampleid}_busco_report.txt" }

    script:
    """
    busco -i $genome \
        -o output-busco-${sampleid} \
        -l bacteria_odb10 \
        --mode genome \
        --cpu ${params.threads} \
        --offline \
        --download_path ${params.busco_db}
    """
}

process RUN_QUAST {
    input:
    tuple val(sampleid), path(genome)

    output:
    tuple val(sampleid), path("output-quast-${sampleid}/report.txt")

    publishDir "${params.output_dir}/output-quast/", mode: 'copy', saveAs: { "${sampleid}_report.txt" }

    script:
    """
    	quast.py $genome -o output-quast-${sampleid}
    """
}

process RUN_BLAST {
    input:
    tuple val(sampleid), path(genome)

    output:
    tuple val(sampleid), path("${sampleid}-blast-genome-vs-db.tsv")

    publishDir "${params.output_dir}/output-blast/", mode: 'copy'

    script:
    """
    blastn \
        -task megablast \
        -query $genome \
        -db ${params.blast_db} \
        -outfmt '6 qseqid staxids bitscore std' \
        -max_target_seqs 1 \
        -max_hsps 1 \
        -num_threads ${params.threads} \
        -evalue 1e-25 \
        -out ${sampleid}-blast-genome-vs-db.tsv
    """
}

process RUN_BLOBTOOLS {
    input:
    tuple val(sampleid), path(genome), path(bam), path(bam_index), path(blast_results)

    output:
    val(sampleid)
    path("${sampleid}_table.tsv")
    path("${sampleid}_blobplot.png")
    path("${sampleid}_blobplot_read_cov.png")

    publishDir "${params.output_dir}/output-blobtools/", mode: 'copy'


    script:
    """

    # prepare blast results
    blobtools taxify -f $blast_results -m ${params.blob_lookup} -s 0 -t 1 -o blast_taxified

    # create database
    #blobtools create -i $genome -b $bam -t $blast_results -o output-blobtools-${sampleid}
    blobtools create -i $genome -b $bam -t blast_taxified.${sampleid}-blast-genome-vs-db.tsv.taxified.out -o output-blobtools-${sampleid} --db ${params.blob_taxdb}
    # produce results table
    blobtools view -i output-blobtools-${sampleid}.blobDB.json -r all -o taxonomy-${sampleid}
    
    # generate figures
    blobtools plot -i output-blobtools-${sampleid}.blobDB.json -r genus

    # rename files
    mv "taxonomy-${sampleid}.output-blobtools-${sampleid}.blobDB.table.txt" "${sampleid}_table.tsv"
    mv "output-blobtools-${sampleid}.blobDB.json.bestsum.genus.p8.span.100.blobplot.bam0.png" "${sampleid}_blobplot.png"
    mv "output-blobtools-${sampleid}.blobDB.json.bestsum.genus.p8.span.100.blobplot.read_cov.bam0.png" "${sampleid}_blobplot_read_cov.png"
    """
}

process RUN_BAKTA {
    input:
    tuple val(sampleid), path(genome)

    output:
    tuple val(sampleid),  file("output-bakta-${sampleid}/${sampleid}.faa"), file("output-bakta-${sampleid}/${sampleid}.gff3")

    publishDir "${params.output_dir}/output-bakta/", mode: 'copy'

    script:
    """

    echo $genome $sampleid

    bakta ${genome} \
	--output output-bakta-${sampleid} \
	--threads ${params.threads} \
	--db ${params.bakta_db} \
	--prefix ${sampleid}
    """
 }

process FINALIZE_RESULTS {
    input:
    tuple val(sample_id), file(fastqc_results), file(contigs), file(quast_results), file(busco_results),

    output:
    file("${sample_id}_final_report.txt") into final_output

    publishDir "./output", mode: 'copy'

    script:
    """
    echo "Sample ID: $sample_id" > ${sample_id}_final_report.txt
    echo "FASTQC Results: $fastqc_results" >> ${sample_id}_final_report.txt
    echo "Contigs: $contigs" >> ${sample_id}_final_report.txt
    echo "QUAST Results: $quast_results" >> ${sample_id}_final_report.txt
    """
}


// Process to Create a Final Summary Table
// process createFinalTable {
//     input:
//     file(final_reports) from final_output.collect()
//
//     output:
//     file("final_results_table.txt")
//
//     script:
//     """
//     # Create the header for the table
//     echo -e "SampleID\tFASTQCResults\tContigs\tQUASTResults" > final_results_table.txt
//
//     # Loop through all final report files and append the necessary data to the table
//     for report in ${final_reports}; do
//         sample_id=\$(basename \$report _final_report.txt)
//         fastqc_result=\$(grep -oP 'FASTQC Results: \K.+' \$report)
//         contigs=\$(grep -oP 'Contigs: \K.+' \$report)
//         quast_result=\$(grep -oP 'QUAST Results: \K.+' \$report)
//
//         echo -e "\$sample_id\t\$fastqc_result\t\$contigs\t\$quast_result" >> final_results_table.txt
//     done
//     """
// }


workflow {
    fastq_files = Channel
    .fromFilePairs(params.reads, flat: false)
    .map { sample_id, files -> tuple(sample_id, files[0], files[1]) }


    // Step 0: Count reads
    read_info = fastq_files | ASSESS_READS
    read_info.view()

    // Step 1: Run FASTQC
    fastqc_report = fastq_files | RUN_FASTQC

    // Step 2: Trim reads with FASTP
    trimmed_fastq = fastq_files | RUN_FASTP

    // Step 3: Assemble genome with SPADES
    genome = trimmed_fastq | RUN_SPADES

    // Step 4.1: Run BWA and Samtools
    bwa_channel = trimmed_fastq.join(genome)
    bam = RUN_BWA(bwa_channel)
    //bam = RUN_BWA(trimmed_fastq, genome)

    // Step 4.2: Run BUSCO
    busco_results = genome | RUN_BUSCO

    // Step 4.3: Run QUAST
    quast_report = genome | RUN_QUAST

    // Step 4.4: Run BLAST
    blast_results = genome | RUN_BLAST

    // Step 5: Run Blobtools
    busco_channel = genome.join(bam).join(blast_results)
    //blob_output = RUN_BLOBTOOLS(genome, bam, blast_results)
    blob_output = RUN_BLOBTOOLS(busco_channel)

//     // Step 6: Contamination checking
//     contamination_check(blob_output)
//
    // Step 7: Annotate genome with BAKTA
    bakta_output = genome | RUN_BAKTA

//     // Final output: Prokka annotated proteome FAA file
//     annotated_proteome.view { it -> println("Final FAA file: ${it}") }

//     fastqc_report.view { println("FASTQC results: $it") }
//     trimmed_fastq.view { println("FASTP results: $it") }
//     genome.view { println("Assembly results: $it") }
//     quast_report.view { println("QUAST report: $it") }

    // Merge results for finalization
//     finalize_input = fastqc_results
//     .zip(assembly) { fastqc, asm -> tuple(fastqc[0], fastqc[1], asm[1]) }
//     .zip(quast_report) { data, quast -> tuple(data[0], data[1], data[2], quast[1]) }


//     finalize_input.view { println "Final input to finalizeResults: $it" }


//     finalize_input | finalizeResults
}

// workflow.onComplete {
//     println "Workflow completed successfully!"
//     println "Reports generated: report.html, timeline.html, trace.txt"
//     println "Final results available in the 'output/' directory."
// }

