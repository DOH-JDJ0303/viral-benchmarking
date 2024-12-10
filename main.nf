#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
=============================================================================================================================
    PIPELINE PARAMETERS
=============================================================================================================================
*/
// Define parameters
params.input            = null
params.outdir           = 'results'
params.length_threshold = 1.2
params.read_limit       = 1000000

def result_file = file(params.outdir).resolve('metrics.csv')

// Check parameters
if( ! params.input ){ exit 1, "Please supply a samplesheet using the '--input' paramater." }
if( ! file(params.input).exists() ){ exit 1, "${params.input} does not exist!" }
if( params.length_threshold < 1 ){ exit 1, "'--length_threshold' must be greater than 1" }

/*
=============================================================================================================================
    MAIN WORKFLOW
=============================================================================================================================
*/
workflow {

    /*
    =============================================================================================================================
        LOAD SAMPLESHEET
    =============================================================================================================================
    */ 
    // Create channel from samplesheet
    Channel.fromPath(params.input)
        .splitCsv(header: true)
        .map{ checkSampleSheet(it) } // Check samplesheet
        .set{ ch_input }
            
    /*
    =============================================================================================================================
        IDENTIFY & EXTRACT TRUTH REGION FROM SAMPLE ASSEMBLY ( IF NEEDED )
        - Determine sequence length of sample and truth assemblies
        - Determine if sample length exceeds the truth length by more than the proportion specified by 'length_threshold'
        - Extract the truth region from the sample using Minimap2 and BedTools. 
    =============================================================================================================================
    */
    // Compare sequence lengths and compare to length threshold
    ch_input
        .map{ [ it.sample, it.assembly, it.truth ] }
        .splitFasta(record: [id: true, seqString: true])
        .splitFasta(record: [id: true, seqString: true])
        .groupTuple(by: 0)
        .map{ sample, sfa, tfa -> [ sample, [ slen: faLength(sfa), tlen: faLength(tfa) ] ] }
        .join( ch_input.map{ [ it.sample, it ] }, by: 0 )
        .map{ sample, lengths, inputs -> inputs + lengths }
        .branch{ it ->
                 normal: it.slen <= params.length_threshold * it.tlen 
                 toolong: it.slen > params.length_threshold * it.tlen 
                }
        .set{ ch_input }

    // Determine overlap between the truth and sample and select the sample contig with the longest alignment length
    // Uses contig name and start/stop location from PAF file to build a BED file for extraction
    MINIMAP2_ASM(
        ch_input.toolong.map{ it -> [ it.sample, it.assembly, it.truth ] }
    )
    MINIMAP2_ASM
        .out
        .paf
        .splitCsv(header: false, sep: '\t')
        .groupTuple(by: 0)
        .map{ sample, paf -> [ sample, paf.find{ line -> line[10] == paf.collect{ it -> it[10] }.max() }[0,2,3] ] }
        .join( ch_input.toolong.map{ [ it.sample, it.assembly ] } )
        .set{ ch_bed }

    // Extract truth region from the sample using BedTools
    // Extracted sequences are combined with samples that did not need extraction (i.e., the 'normal' samples)
    BEDTOOLS(
        ch_bed
    )
    BEDTOOLS
        .out
        .fa
        .join( ch_input.toolong.map{ it -> [ it.sample, it ] } )
        .map{ sample, fragment, flen, it -> it + [ fragment: fragment, flen: flen ] }
        .concat( ch_input.normal.map{ it + [ fragment: null, flen: null ] } )
        .set{ ch_input }

    /*
    =============================================================================================================================
        COMPARE TRUTH AND SAMPLE ASSEMBLIES
        - Gather initial metrics using Nextclade
        - Capture any Nextclade errors
        - Account for unaligned sequences on the terminal ends (Nextlade ignores these)
        - Calculate accuracy and completeness
    =============================================================================================================================
    */
    // Gather initial metrics with Nextclade
    NEXTCLADE(
        ch_input.map{ [ it.sample, it.fragment ? it.fragment : it.assembly, it.truth ] }
    )
    // Load Nextclade results
    NEXTCLADE
        .out
        .json
        .splitJson()
        .set{ ch_nc }
    // Capture Nextclade errors
    ch_nc
        .filter{ sample, item -> item.key == 'errors' }
        .map{ sample, item -> [ sample, item.value[0] ? item.value[0].errors[0] : null ]  }
        .set{ ch_nc_errors }
    // Extract Nextclade results and gather assembly metrics
    ch_nc
        .filter{ sample, item -> item.key == 'results' }
        .join( ch_input.map{ [ it.sample, it.slen, it.tlen ] }, by: 0 )
        .join( ch_nc_errors, by: 0 )
        .map{ sample, results, slen, tlen, nc_errors -> [ sample: sample, nc_errors: nc_errors, slen: slen ] + gatherMetrics(results + [ tlen: tlen ]) }
        .set{ ch_assembly_metrics }

    /*
    =============================================================================================================================
        GATHER READ-LEVEL METRICS (OPTIONAL)
    =============================================================================================================================
    */
    // Assign read channel
    ch_input
        .filter{ it.fastq_1 }
        .map{ [ sample: it.sample, reads: [ it.fastq_1, it.fastq_2 ? it.fastq_2 : null ], read_type: 'short' ] }
        .concat( ch_input.filter{ it.fastq_l }.map{ [ sample: it.sample, reads: it.fastq_l, read_type: 'long' ] } )
        .set{ ch_reads }
    // Down-sample reads, if needed
    ch_reads
        .map{ [ it, it.read_type == 'short' ? it.reads[0] : it.reads ] }
        .splitFastq(record: true)
        .groupTuple(by: 0)
        .map{ it, fq -> it + [ n_reads: it.read_type == "short" ? fq.size()*2 : fq.size() ] }
        .branch{ it ->
            ok: it.n_reads <= params.read_limit
            toomuch: it.n_reads > params.read_limit }
        .set{ ch_reads }
    // Run 'seqtk sample' for down-sampling
    SEQTK_SAMPLE( 
        ch_reads.toomuch.map{ [it.sample, it.reads, it.read_type ] }
    )
    // Combine down-sampled channel with non-down-sampled channel
    SEQTK_SAMPLE
        .out
        .reads
        .map{ sample, reads, read_type -> [ sample: sample, reads: reads, read_type: read_type, n_reads: params.read_limit ] }
        .concat(ch_reads.ok)
        .map{ [ it.sample, it ] }
        .join( ch_input.map{ [ it.sample, [ truth: it.truth ] ] }, by: 0 )
        .join( ch_reads.ok.concat( ch_reads.toomuch ).map{ [ it.sample, [ n_reads_full: it.n_reads ] ] }, by: 0 )
        .map{ sample, it, truth, n_reads_full -> it + truth + n_reads_full }
        .map{ it + [ read_sample_rate: it.n_reads / it.n_reads_full ] }
        .set{ ch_reads }
    // Align short reads to the truth assembly with 'bwa mem'
    BWA_MEM(
        ch_reads.filter{ it.read_type == 'short' }.map{ [ it.sample, it.reads, it.truth ] }
    )
    // Align long reeads to the truth assembly with 'minimap2'
    MINIMAP2_MAP(
        ch_reads.filter{ it.read_type == 'long' }.map{ [ it.sample, it.reads, it.truth ] }
    )
    // Combine the short and long mapped reads into one channel
    BWA_MEM
        .out
        .aln
        .concat(MINIMAP2_MAP.out.aln)
        .set{ ch_aln }
    // Get read alignment stats with 'samtools coverage'
    SAMTOOLS_COVERAGE(
        ch_aln
    )
    // Gather read-level metrics from the samtools output
    SAMTOOLS_COVERAGE
        .out
        .coverage
        .splitCsv(header: true, sep: '\t')
        .join( ch_reads.map{ [ it.sample, it ] }, by: 0 )
        .map{ sample, coverage, read_info -> read_info + coverage }
        .map{ [ sample: it.sample, 
                n_reads: it.n_reads_full, 
                per_reads_mapped: ( 100 * it.numreads.toFloat() / it.n_reads ).round(),
                read_cov: it.coverage.toFloat(),
                mean_read_depth: ( it.meandepth.toFloat() / it.read_sample_rate ).round(), 
                mean_read_q: it.meanbaseq.toFloat(), 
                mean_map_q: it.meanmapq.toFloat() ] }
        .set{ ch_read_metrics }

    /*
    =============================================================================================================================
        MERGE METRIC CHANNELS
    =============================================================================================================================
    */
    // Merge the assembly and read metric channels
    ch_assembly_metrics
        .map{ [ it.sample, it ] }
        .join( ch_read_metrics.map{ [ it.sample, it ] }, by: 0, remainder: true )
        .map{ sample, assembly_metrics, read_metrics -> assembly_metrics + ( read_metrics instanceof LinkedHashMap ? read_metrics : [ ] ) }
        .set{ ch_all_metrics }
    // Get a list of linked hashmap keys across all samples and fill any missing keys in each sample with 'null'
    ch_all_metrics
        .map{ it.keySet() }
        .collect()
        .map{ [ it.flatten().unique() ] }
        .combine( ch_all_metrics )
        .map{ cols, metrics -> fillMissingKeys(cols, metrics) }
        .set{ ch_all_metrics }
    // Export the metrics to a file
    ch_all_metrics
        .take(1)
        .map{ it.keySet().join(',') }
        .concat(ch_all_metrics.map{ it.values().join(',') })
        .collectFile(name: 'result.csv', sort: 'index', newLine: true)
        .subscribe{ file(it).copyTo(result_file) }
}

// Actions upon workflow completion
workflow.onComplete {
    println "\nWorkflow complete!\nResults can be found at ${result_file}"
}

/*
=============================================================================================================================
    FUNCTIONS
=============================================================================================================================
*/
// Function for checking and formatting the input samplesheet
def checkSampleSheet(LinkedHashMap row){
    // Check for missing columns
    def missing_cols = [ "sample", "assembly", "truth" ].findAll { !row.containsKey(it) }
    if( missing_cols ){ exit 1, "The supplied samplesheet does not contain one or more required column: ${missing_cols.join(', ')}" }
    
    // Check if input files exist
    def assembly = file(row.assembly, checkIfExists: true)
    def truth = file(row.truth, checkIfExists: true)
    def fastq_1 = row.fastq_1 ? file(row.fastq_1, checkIfExists: true) : null
    def fastq_2 = row.fastq_2 ? file(row.fastq_2, checkIfExists: true) : null
    def fastq_l = row.fastq_l ? file(row.fastq_l, checkIfExists: true) : null

    return [ sample: row.sample, assembly: assembly, truth: truth, fastq_1: fastq_1, fastq_2: fastq_2, fastq_l: fastq_l ]
}

// Function for getting the length of all contigs in a fasta file
def faLength(contigs){
    return contigs.unique().collect{ it.seqString.length() }.sum()
}

// Function for gathering assembly metrics
def gatherMetrics(LinkedHashMap row){
    def metrics            = row.value[0] ? row.value[0] : null
    def totalSubstitutions = metrics ? metrics.totalSubstitutions : null
    def totalDeletions     = metrics ? metrics.totalDeletions : null
    def totalInsertions    = metrics ? metrics.totalInsertions : null
    def totalMissing       = metrics ? metrics.totalMissing : null
    def totalNonACGTNs     = metrics ? metrics.totalNonACGTNs : null
    def terminalMiss       = metrics ? ( metrics.nucleotideComposition.'-' ? metrics.nucleotideComposition.'-' : 0 ) : null

    def alignLen           = metrics ? row.tlen - terminalMiss - totalMissing - totalNonACGTNs : null
    def agreement          = metrics ? alignLen - totalSubstitutions - totalDeletions - totalInsertions : null
    def completeness       = metrics ? 100 * alignLen / row.tlen : null
    def accuracy           = metrics ? 100 * agreement / alignLen : null

    return [ tlen: row.tlen,
             sub: totalSubstitutions, 
             del: totalDeletions, 
             ins: totalInsertions, 
             miss: totalMissing, 
             other: totalNonACGTNs, 
             term: terminalMiss, 
             completeness: completeness,
             accuracy: accuracy ]
}

// Function for filling in missing linked hashmap keys with null
def fillMissingKeys(cols, row){
    def missing_keys = cols.findAll { !row.containsKey(it) }
    if(missing_keys.size() > 0){
        def makeup_keys = [:]
        missing_keys.each{ makeup_keys[it] = null }
        row = row + makeup_keys
    }

    return row
}

/*
=============================================================================================================================
    PROCESSES
=============================================================================================================================
*/

process MINIMAP2_ASM {
    input:
    tuple val(sample), path(assembly), path(truth)

    output:
    tuple val(sample), path("*.paf"),   emit: paf

    shell:
    """
    minimap2 -x asm20 ${truth} ${assembly} > ${sample}.paf
    """
}

process BEDTOOLS {
    input:
    tuple val(sample), val(bed), path(assembly)

    output:
    tuple val(sample), path("*.fa"), val("${bed[2]-bed[0]}"), emit: fa

    script:
    """
    echo "${bed.join('\t')}" > coords.bed
    bedtools getfasta -fi ${assembly} -fo "${sample}.${bed[1]}-${bed[2]}.fa" -bed coords.bed
    """
}

process NEXTCLADE {
    input:
    tuple val(sample), path(assembly), path(truth)

    output:
    tuple val(sample), path("results/*.fasta"),   emit: fasta
    tuple val(sample), path("results/*.json"), emit: json
    
    script:
    """
    nextclade run -r ${truth} -O results ${assembly}
    """
}

process SEQTK_SAMPLE {
    input:
    tuple val(sample), path(reads), val(read_type)

    output:
    tuple val(sample), path("*.fastq.gz"), val(read_type), emit: reads

    

    script:
    if( read_type == 'short' ){
        cmd = "seqtk sample ${reads[0]} ${params.read_limit} > ${sample}_R1.fastq"
        if( { reads[1] } ){
            cmd = "${cmd} && seqtk sample ${reads[0]} ${params.read_limit} > ${sample}_R2.fastq"
        }
    }else{
        cmd = "seqtk sample ${reads} ${params.read_limit} > ${sample}_long.fastq"
    }
    """
    $cmd
    gzip *.fastq
    """

}

process BWA_MEM {
    input:
    tuple val(sample), path(reads), path(truth)

    output:
    tuple val(sample), path('*bam'), emit: aln

    script:
    """
    # setup for pipe
    set -euxo pipefail

    # index the reference
    bwa index ${truth}

    # run bwa mem, select only mapped reads, convert to .bam, and sort
    bwa mem -t ${task.cpus} ${truth} ${reads[0]} ${reads[1] ? reads[1] : ''} | samtools view -b -F 4 - | samtools sort - > ${sample}.bam
    """
}

process MINIMAP2_MAP {
    input:
    tuple val(sample), path(reads), path(truth)

    output:
    tuple val(sample), path('*bam'), emit: aln

    script:
    """
    # setup for pipe
    set -euxo pipefail

    # run bwa mem, select only mapped reads, convert to .bam, and sort
    minimap2 -a -x map-ont -t ${task.cpus} ${truth} ${reads} | samtools view -b -F 4 - | samtools sort - > ${sample}.bam
    """
}

process SAMTOOLS_COVERAGE {
    input:
    tuple val(sample), path(bam)

    output:
    tuple val(sample), path("*.coverage.txt"), emit: coverage

    script:
    """ 
    # gather read stats
    samtools coverage ${sample}.bam > ${sample}.coverage.txt
    """
}