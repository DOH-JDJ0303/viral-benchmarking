#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
=============================================================================================================================
    PIPELINE PARAMETERS
=============================================================================================================================
*/
// Define results file
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
        .map{ sample, sfa, tfa -> [ sample, [ sampleLength: faLength(sfa), truthLength: faLength(tfa) ] ] }
        .join( ch_input.map{ [ it.sample, it ] }, by: 0 )
        .map{ sample, lengths, inputs -> inputs + lengths }
        .branch{ it ->
                 normal: it.sampleLength <= params.length_threshold * it.truthLength 
                 toolong: it.sampleLength > params.length_threshold * it.truthLength 
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
        - Account for "missing bases" introduced by the reference genome
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
        .set{ ch_ncErrors }
    // Extract Nextclade results and gather assembly metrics
    ch_nc
        .filter{ sample, item -> item.key == 'results' }
        .join( ch_input.map{ [ it.sample, it.sampleLength, it.truthLength ] }, by: 0 )
        .join( ch_ncErrors, by: 0 )
        .map{ sample, results, sampleLength, truthLength, ncErrors -> [ sample: sample, ncErrors: ncErrors, sampleLength: sampleLength ] + gatherMetrics(results + [ truthLength: truthLength ]) }
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
                nReads: it.n_reads_full, 
                perReadsMapped: ( 100 * it.numreads.toFloat() / it.n_reads ).round(),
                readCov: it.coverage.toFloat(),
                meanReadDepth: ( it.meandepth.toFloat() / it.read_sample_rate ).round(), 
                meanReadQual: it.meanbaseq.toFloat(), 
                meanMapQual: it.meanmapq.toFloat() ] }
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
        .map{ it.assembly = it.sample.seq
              it.sample = it.sample.id
              it 
        }
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
    def truth    = file(row.truth, checkIfExists: true)
    def fastq_1  = row.fastq_1 ? file(row.fastq_1, checkIfExists: true) : null
    def fastq_2  = row.fastq_2 ? file(row.fastq_2, checkIfExists: true) : null
    def fastq_l  = row.fastq_l ? file(row.fastq_l, checkIfExists: true) : null
    def sample   = [ id: row.sample, seq: assembly.getSimpleName() ]

    return [ sample: sample, assembly: assembly, truth: truth, fastq_1: fastq_1, fastq_2: fastq_2, fastq_l: fastq_l ]
}

// Function for getting the length of all contigs in a fasta file
def faLength(contigs){
    return contigs.unique().collect{ it.seqString.length() }.sum()
}

// Function for gathering assembly metrics
def gatherMetrics(LinkedHashMap row){
    // Define metrics
    def metrics            = row.value[0] ? row.value[0] : null

    /*
    ==============================================================
        Default Nextclade Metrics
    ==============================================================
    */        
    def totalSubstitutions = metrics ? metrics.totalSubstitutions : null
    def totalDeletions     = metrics ? metrics.totalDeletions : null
    def totalNonACGTNs     = metrics ? metrics.totalNonACGTNs : null

    /*
    ==============================================================
        Modified Nextclade Metrics
    ==============================================================
    */    
    // Missing data - Nextclade does not count sites in the sample that are missing from the truth termini
    def termMissing        = metrics ? ( metrics.nucleotideComposition.'-' ? metrics.nucleotideComposition.'-' : 0 ) : null
    def totalMissing       = metrics ? metrics.totalMissing + termMissing : null
    // False insertions due to differences between the truth sequence and the reference used to create the assembly
    def fauxInsertions     = metrics ? 0 : null
    if( metrics ){ metrics.insertions.each{ fauxInsertions = fauxInsertions + it.ins.count( 'N' ) + it.ins.count( '-' ) } }
    def totalInsertions    = metrics ? metrics.totalInsertions - fauxInsertions : null

    /*
    ==============================================================
        New metrics
    ==============================================================
    */
    // Alignment length - number of comparable sites between the truth and sample
    def alignLen           = metrics ? row.truthLength - totalMissing : null
    // Completeness     - percentage of the truth sequence captured by the sample
    def completeness       = metrics ? 100 * alignLen / row.truthLength : null
    // Accuracy         - percentage of agreements between the sample and the truth at comparable sites
    def disagreement       = metrics ? totalSubstitutions + totalDeletions + totalInsertions + totalNonACGTNs : null
    def agreement          = metrics ? alignLen - disagreement : null
    def accuracy           = metrics ? 100 * agreement / alignLen : null

    /*
    ==============================================================
        Final Metrics
    ==============================================================
    */    
    return [ truthLength: row.truthLength,
             substitutions: totalSubstitutions,
             deletions: totalDeletions, 
             insertions: totalInsertions,
             fauxInsertions: fauxInsertions,
             missing: totalMissing,
             terminiMissing: termMissing,
             nonACGTNs: totalNonACGTNs,
             aligned: alignLen,
             disagreements: disagreement,
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
    container 'docker.io/staphb/minimap2:2.28'

    input:
    tuple val(sample), path(assembly), path(truth)

    output:
    tuple val(sample), path("*.paf"),   emit: paf

    script:
    prefix = "${sample.id}_${sample.seq}"
    """
    minimap2 -x asm20 ${truth} ${assembly} > ${prefix}.paf
    """
}

process BEDTOOLS {
    container 'docker.io/staphb/bedtools:2.31.1'

    input:
    tuple val(sample), val(bed), path(assembly)

    output:
    tuple val(sample), path("*.fa"), val("${bed[2]-bed[0]}"), emit: fa

    script:
    prefix = "${sample.id}_${sample.seq}"
    """
    echo "${bed.join('\t')}" > coords.bed
    bedtools getfasta -fi ${assembly} -fo "${prefix}.${bed[1]}-${bed[2]}.fa" -bed coords.bed
    """
}

process NEXTCLADE {
    container 'docker.io/nextstrain/nextclade:3.9.1'
    publishDir "${params.outdir}/nextclade/${prefix}/", mode: 'symlink'

    input:
    tuple val(sample), path(assembly), path(truth)

    output:
    tuple val(sample), path("results/*.fasta"), emit: fasta
    tuple val(sample), path("results/*.json"),  emit: json
    
    script:
    prefix = "${sample.id}_${sample.seq}"
    """
    nextclade run -r ${truth} -O results ${assembly}
    """
}

process SEQTK_SAMPLE {
    container 'docker.io/staphb/seqtk:1.4'

    input:
    tuple val(sample), path(reads), val(read_type)

    output:
    tuple val(sample), path("*.fastq.gz"), val(read_type), emit: reads

    script:
    if( read_type == 'short' ){
        cmd = "seqtk sample ${reads[0]} ${params.read_limit} > ${sample.id}_R1.fastq"
        if( { reads[1] } ){
            cmd = "${cmd} && seqtk sample ${reads[0]} ${params.read_limit} > ${sample.id}_R2.fastq"
        }
    }else{
        cmd = "seqtk sample ${reads} ${params.read_limit} > ${sample.id}_long.fastq"
    }
    """
    $cmd
    gzip *.fastq
    """

}

process BWA_MEM {
    container 'quay.io/biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:219b6c272b25e7e642ae3ff0bf0c5c81a5135ab4-0'

    input:
    tuple val(sample), path(reads), path(truth)

    output:
    tuple val(sample), path('*bam'), emit: aln

    script:
    prefix = "${sample.id}_${sample.seq}"
    """
    # setup for pipe
    set -euxo pipefail

    # index the reference
    bwa index ${truth}

    # run bwa mem, select only mapped reads, convert to .bam, and sort
    bwa mem -t ${task.cpus} ${truth} ${reads[0]} ${reads[1] ? reads[1] : ''} | samtools view -b -F 4 - | samtools sort - > ${prefix}.bam
    """
}

process MINIMAP2_MAP {
    container 'docker.io/staphb/minimap2:2.28'

    input:
    tuple val(sample), path(reads), path(truth)

    output:
    tuple val(sample), path('*bam'), emit: aln

    script:
    prefix = "${sample.id}_${sample.seq}"
    """
    # setup for pipe
    set -euxo pipefail

    # run bwa mem, select only mapped reads, convert to .bam, and sort
    minimap2 -a -x map-ont -t ${task.cpus} ${truth} ${reads} | samtools view -b -F 4 - | samtools sort - > ${prefix}.bam
    """
}

process SAMTOOLS_COVERAGE {
    container 'docker.io/staphb/samtools:1.21'

    input:
    tuple val(sample), path(bam)

    output:
    tuple val(sample), path("*.coverage.txt"), emit: coverage

    script:
    prefix = "${sample.id}_${sample.seq}"
    """ 
    # gather read stats
    samtools coverage ${bam} > ${prefix}.coverage.txt
    """
}