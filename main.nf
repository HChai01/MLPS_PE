nextflow.enable.dsl = 2

// ========================================================================================
// Workflow Imports
// ========================================================================================

include { MLPS_PE_PIPELINE } from './subworkflows/mlps_pe'


// ========================================================================================
// Parameter Parsing Helpers
// ========================================================================================

// Parse yes/no style CLI flags into strict booleans before the workflow starts.
def parseBoolParam(value, name) {
    if( value instanceof Boolean ) {
        return value
    }
    def normalized = value?.toString()?.trim()?.toLowerCase()
    if( normalized in ['yes', 'true', '1'] ) {
        return true
    }
    if( normalized in ['no', 'false', '0'] ) {
        return false
    }
    throw new IllegalArgumentException("Parameter --${name} must be yes/no or true/false. Received: ${value}")
}

// Normalize human-readable byte sizes so the Python cleanup helper receives a concrete integer threshold regardless of how the user expressed the value.
def parseSizeParam(value, name) {
    def normalized = value?.toString()?.trim()?.toUpperCase()
    if( !normalized ) {
        throw new IllegalArgumentException("Parameter --${name} must not be empty")
    }

    def matcher = normalized =~ /^(\d+)(B|KB|MB|GB)?$/
    if( !matcher.matches() ) {
        throw new IllegalArgumentException("Parameter --${name} must look like 1024, 1KB, 1MB, or 1GB. Received: ${value}")
    }

    long amount = matcher[0][1] as Long
    def unit = matcher[0][2] ?: 'B'
    def multiplier = [
        'B' : 1L,
        'KB': 1024L,
        'MB': 1024L * 1024L,
        'GB': 1024L * 1024L * 1024L,
    ][unit]

    return amount * multiplier
}

// Canonicalize path-style parameters up front so the workflow logs and all downstream helpers see the same resolved filesystem location.
def resolvePathParam(value, name, boolean mustExist = true) {
    if( value == null ) {
        if( mustExist ) {
            throw new IllegalArgumentException("Missing required path parameter: --${name}")
        }
        return null
    }

    def resolved = new File(value.toString()).getCanonicalFile()
    if( mustExist && !resolved.exists() ) {
        throw new IllegalArgumentException("Path provided via --${name} was not found: ${value}")
    }
    return resolved.absolutePath
}

def resolveProjectNameFromOutdir(outdirValue) {
    if( outdirValue == null ) {
        return null
    }

    def resolvedOutdir = new File(outdirValue.toString()).getCanonicalFile()
    return resolvedOutdir.name
}

def isValidPrimerFasta(value) {
    if( value == null ) {
        return false
    }

    def normalized = value.toString().trim()
    if( !normalized || normalized.equalsIgnoreCase('null') ) {
        return false
    }

    def primerFile = new File(normalized).getCanonicalFile()
    if( !primerFile.exists() || !primerFile.isFile() ) {
        return false
    }

    boolean hasHeader = false
    boolean hasSequence = false
    primerFile.eachLine { rawLine ->
        if( hasSequence ) {
            return
        }

        def line = rawLine.trim()
        if( !line ) {
            return
        }

        if( !hasHeader ) {
            hasHeader = line.startsWith('>')
        }
        else if( !line.startsWith('>') ) {
            hasSequence = true
        }
    }

    return hasHeader && hasSequence
}

// Resolve primer files for different input
def resolvePrimerFastaParams(fwdValue, revValue) {
    def dummyPrimer = new File(projectDir.toString(), 'reference/dummy_primer.fasta').getCanonicalFile()
    if( !dummyPrimer.exists() || !dummyPrimer.isFile() ) {
        throw new IllegalArgumentException("Default dummy primer FASTA was not found: ${dummyPrimer.absolutePath}")
    }

    def resolvedFwd = fwdValue ? new File(fwdValue.toString()).getCanonicalFile().absolutePath : null
    def resolvedRev = revValue ? new File(revValue.toString()).getCanonicalFile().absolutePath : null
    def fwdValid = isValidPrimerFasta(fwdValue)
    def revValid = isValidPrimerFasta(revValue)

    if( fwdValid && revValid ) {
        return [fwd_primer_seq: resolvedFwd, rev_primer_seq: resolvedRev]
    }

    if( fwdValid && !revValid ) {
        log.warn "Parameter --rev_primer_seq is missing or invalid. Falling back to --fwd_primer_seq for both primer inputs."
        return [fwd_primer_seq: resolvedFwd, rev_primer_seq: resolvedFwd]
    }

    if( !fwdValid && revValid ) {
        log.warn "Parameter --fwd_primer_seq is missing or invalid while --rev_primer_seq is valid. Falling back to ${dummyPrimer.absolutePath} for both primer inputs."
        return [fwd_primer_seq: dummyPrimer.absolutePath, rev_primer_seq: dummyPrimer.absolutePath]
    }

    log.warn "Primer FASTA inputs are missing or invalid. Falling back to ${dummyPrimer.absolutePath} for both primer inputs."
    return [fwd_primer_seq: dummyPrimer.absolutePath, rev_primer_seq: dummyPrimer.absolutePath]
}


// ========================================================================================
// Workflow Finalizers
// ========================================================================================

// Always attempt to merge the published summary TSVs into one project-level table when the workflow terminates, even if execution ended in failure after some samples passed.
def mergePublishedSummaryResults() {
    def outdirValue = params.outdir?.toString()
    def projectNameValue = resolveProjectNameFromOutdir(outdirValue)
    if( !outdirValue || !projectNameValue ) {
        log.warn "Skipping final summary merge because outdir/project_name was not resolved"
        return
    }

    def summaryDir = new File(outdirValue, 'summary')
    if( !summaryDir.exists() ) {
        log.warn "Skipping final summary merge because summary directory does not exist yet: ${summaryDir.absolutePath}"
        return
    }

    def mergeScript = new File(projectDir.toString(), 'python/merge_tsv.py')
    def inputCsvPath = new File(params.inputcsv.toString()).getCanonicalPath()
    def outputTsv = new File(summaryDir, "${projectNameValue}.results.tsv")
    def command = [
        'python3',
        mergeScript.absolutePath,
        '--indir', summaryDir.absolutePath,
        '--id', 'job_id',
        '--samplecsv', inputCsvPath,
        '--outputtsv', outputTsv.absolutePath
    ]

    log.info "Running final summary merge: ${command.join(' ')}"

    def process = new ProcessBuilder(command)
        .directory(new File(projectDir.toString()))
        .redirectErrorStream(true)
        .start()

    def outputLines = []
    process.inputStream.withReader { reader ->
        reader.eachLine { line -> outputLines << line }
    }

    def exitCode = process.waitFor()
    if( outputLines ) {
        outputLines.each { line -> log.info "[merge_tsv] ${line}" }
    }

    if( exitCode != 0 ) {
        log.warn "Final summary merge exited with code ${exitCode}; partial summary TSVs remain in ${summaryDir.absolutePath}"
        return
    }

    if( outputTsv.exists() ) {
        log.info "Merged project summary written to ${outputTsv.absolutePath}"
    }
    else {
        log.warn "Final summary merge completed but did not create ${outputTsv.absolutePath}"
    }
}

// Export the executed Nextflow task scripts into a project-level Markdown audit file.
def exportPipelineCommandAudit() {
    def outdirValue = params.outdir?.toString()
    def projectNameValue = resolveProjectNameFromOutdir(outdirValue)
    if( !outdirValue || !projectNameValue ) {
        log.warn "Skipping command audit export because outdir/project_name was not resolved"
        return
    }

    def pipelineInfoDir = new File(outdirValue, 'pipeline_info')
    pipelineInfoDir.mkdirs()

    def traceFile = new File(pipelineInfoDir, 'trace.txt')
    def workDirPath = workflow.workDir?.toString() ?: new File(projectDir.toString(), 'work').absolutePath
    def workDirFile = new File(workDirPath)
    def exportScript = new File(projectDir.toString(), 'python/export_command_md.py')
    def outputMd = new File(pipelineInfoDir, "${projectNameValue}.command.md")

    if( !workDirFile.exists() ) {
        log.warn "Skipping command audit export because work directory does not exist: ${workDirFile.absolutePath}"
        return
    }

    def command = [
        'python3',
        exportScript.absolutePath,
        '--trace', traceFile.absolutePath,
        '--workdir', workDirFile.absolutePath,
        '--outputmd', outputMd.absolutePath,
        '--project', projectNameValue
    ]

    log.info "Running command audit export: ${command.join(' ')}"

    def process = new ProcessBuilder(command)
        .directory(new File(projectDir.toString()))
        .redirectErrorStream(true)
        .start()

    def outputLines = []
    process.inputStream.withReader { reader ->
        reader.eachLine { line -> outputLines << line }
    }

    def exitCode = process.waitFor()
    if( outputLines ) {
        outputLines.each { line -> log.info "[command_audit] ${line}" }
    }

    if( exitCode != 0 ) {
        log.warn "Command audit export exited with code ${exitCode}"
        return
    }

    if( outputMd.exists() ) {
        log.info "Pipeline command audit written to ${outputMd.absolutePath}"
    }
    else {
        log.warn "Command audit export completed but did not create ${outputMd.absolutePath}"
    }
}

// Truncate oversized regular files inside traced Nextflow work directories.
// This only touches files under the real work directory and skips symlinks.
def cleanupTemporaryTaskData() {
    if( !params.cleantmpdata ) {
        log.info "Skipping temporary task data cleanup because --cleantmpdata is disabled"
        return
    }

    def outdirValue = params.outdir?.toString()
    if( !outdirValue ) {
        log.warn "Skipping temporary task data cleanup because outdir was not resolved"
        return
    }

    def pipelineInfoDir = new File(outdirValue, 'pipeline_info')
    def traceFile = new File(pipelineInfoDir, 'trace.txt')
    def workDirPath = workflow.workDir?.toString() ?: new File(projectDir.toString(), 'work').absolutePath
    def workDirFile = new File(workDirPath)
    def cleanupScript = new File(projectDir.toString(), 'python/cleanup_tmp_data.py')

    if( !traceFile.exists() ) {
        log.warn "Skipping temporary task data cleanup because trace file was not found: ${traceFile.absolutePath}"
        return
    }
    if( !workDirFile.exists() ) {
        log.warn "Skipping temporary task data cleanup because work directory does not exist: ${workDirFile.absolutePath}"
        return
    }

    def command = [
        'python3',
        cleanupScript.absolutePath,
        '--trace', traceFile.absolutePath,
        '--workdir', workDirFile.absolutePath,
        '--max-bytes', params.cleantmpdatasize_bytes.toString()
    ]

    log.info "Running temporary task data cleanup for traced work directories with threshold ${params.cleantmpdatasize} (${params.cleantmpdatasize_bytes} bytes)"

    def process = new ProcessBuilder(command)
        .directory(new File(projectDir.toString()))
        .redirectErrorStream(true)
        .start()

    def outputLines = []
    process.inputStream.withReader { reader ->
        reader.eachLine { line -> outputLines << line }
    }

    def exitCode = process.waitFor()
    if( outputLines ) {
        outputLines.each { line -> log.info "[cleanup_tmp] ${line}" }
    }

    if( exitCode != 0 ) {
        log.warn "Temporary task data cleanup exited with code ${exitCode}"
    }
}

// Print a terminal summary directly to stdout.
def printCompletionSummary() {
    def stats = workflow.stats
    def completedAt = workflow.complete ?: new Date()
    def status = workflow.success ? 'SUCCEEDED' : 'FAILED'

    println ''
    println '============================================================'
    println "Pipeline status: ${status}"
    def timeFormat = java.time.format.DateTimeFormatter.ofPattern('dd-MMM-yyyy HH:mm:ss')
    println "Completed at   : ${completedAt.format(timeFormat)}"
    println "Duration       : ${workflow.duration}"
    println "Succeeded      : ${stats.succeededCount}"
    println "Failed         : ${stats.failedCount}"
    println "============================================================"
}


// ========================================================================================
// Input Discovery Helpers
// ========================================================================================

// Parse one CSV line while preserving commas embedded inside quoted fields so the workflow can accept sample sheets exported from spreadsheet applications.
def parseCsvLine(String line) {
    def values = []
    def current = new StringBuilder()
    boolean inQuotes = false
    boolean skipNext = false

    line.toCharArray().eachWithIndex { ch, i ->
        if( skipNext ) {
            skipNext = false
            return
        }

        if( ch == '"' ) {
            if( inQuotes && i + 1 < line.length() && line.charAt(i + 1) == '"' ) {
                current.append('"')
                skipNext = true
            }
            else {
                inQuotes = !inQuotes
            }
            return
        }

        if( ch == ',' && !inQuotes ) {
            values << current.toString().trim()
            current.setLength(0)
            return
        }

        current.append(ch)
    }

    values << current.toString().trim()
    return values
}

def resolveInputReadPath(String value, String columnName, int rowNumber, File baseDir) {
    // Resolve relative FASTQ paths against the sample-sheet directory so the same input CSV works whether the workflow is launched locally or from SBATCH wrappers.
    if( !value ) {
        throw new IllegalArgumentException("Row ${rowNumber} in --inputcsv is missing ${columnName}")
    }

    def candidate = new File(value)
    def resolved = (candidate.isAbsolute() ? candidate : new File(baseDir, value)).getCanonicalFile()
    if( !resolved.exists() || !resolved.isFile() ) {
        throw new IllegalArgumentException("Row ${rowNumber} in --inputcsv references a missing file for ${columnName}: ${value}")
    }
    return resolved
}

def discoverInputPairs(String inputCsvPath) {
    // Ignore blank and comment lines so operational notes can live beside the sample sheet.
    def inputCsv = new File(inputCsvPath)
    def csvDir = inputCsv.getCanonicalFile().parentFile
    def lines = inputCsv.readLines('UTF-8')
        .findAll { line -> line?.trim() && !line.trim().startsWith('#') }

    if( !lines ) {
        throw new IllegalArgumentException("Input CSV is empty: ${inputCsvPath}")
    }

    def header = parseCsvLine(lines[0])
    def requiredHeaders = ['job_id', 'sample_name', 'r1_path', 'r2_path']
    def missingHeaders = requiredHeaders.findAll { !header.contains(it) }
    if( missingHeaders ) {
        throw new IllegalArgumentException(
            "Input CSV is missing required header column(s): ${missingHeaders.join(', ')}"
        )
    }

    def headerIndex = header.withIndex().collectEntries { value, index -> [(value): index] }
    def samples = []

    lines.drop(1).eachWithIndex { line, offset ->
        def rowNumber = offset + 2
        def fields = parseCsvLine(line)
        if( fields.size() != header.size() ) {
            throw new IllegalArgumentException(
                "Row ${rowNumber} in --inputcsv has ${fields.size()} column(s); expected ${header.size()}"
            )
        }

        def jobId = fields[headerIndex.job_id]?.trim()
        def sampleName = fields[headerIndex.sample_name]?.trim()
        if( !jobId ) {
            throw new IllegalArgumentException("Row ${rowNumber} in --inputcsv is missing job_id")
        }
        if( !sampleName ) {
            throw new IllegalArgumentException("Row ${rowNumber} in --inputcsv is missing sample_name")
        }

        def r1File = resolveInputReadPath(fields[headerIndex.r1_path], 'r1_path', rowNumber, csvDir)
        def r2File = resolveInputReadPath(fields[headerIndex.r2_path], 'r2_path', rowNumber, csvDir)

        samples << [
            id               : jobId,
            sample_name      : sampleName,
            r1_original_name : r1File.name,
            r2_original_name : r2File.name,
            r1_canonical_name: "${jobId}.R1.fastq.gz",
            r2_canonical_name: "${jobId}.R2.fastq.gz",
            r1_path          : r1File.absolutePath,
            r2_path          : r2File.absolutePath,
        ]
    }

    if( !samples ) {
        throw new IllegalArgumentException("Input CSV contains no sample rows: ${inputCsvPath}")
    }

    def duplicateIds = samples
        .groupBy { it.id }
        .findAll { sampleId, rows -> rows.size() > 1 }
        .keySet()

    if( duplicateIds ) {
        throw new IllegalArgumentException("Duplicate job_id values detected in --inputcsv: ${duplicateIds.join(', ')}")
    }

    return samples
}


// ========================================================================================
// Workflow Validation
// ========================================================================================

// Validate user parameters early so the workflow fails with a clear message before expensive processing begins.
def validateParams() {
    params.inputcsv = resolvePathParam(params.inputcsv, 'inputcsv')
    params.outdir = resolvePathParam(params.outdir ?: 'results', 'outdir', false)
    params.target_reference = resolvePathParam(params.target_reference, 'target_reference')
    def resolvedPrimers = resolvePrimerFastaParams(params.fwd_primer_seq, params.rev_primer_seq)
    params.fwd_primer_seq = resolvedPrimers.fwd_primer_seq
    params.rev_primer_seq = resolvedPrimers.rev_primer_seq

    params.dedup = parseBoolParam(params.dedup, 'dedup')
    params.remove_human = parseBoolParam(params.remove_human, 'remove_human')
    params.cleantmpdata = parseBoolParam(params.cleantmpdata, 'cleantmpdata')
    params.cleantmpdatasize_bytes = parseSizeParam(params.cleantmpdatasize, 'cleantmpdatasize')

    def genomeType = params.target_genome_type?.toString()?.trim()?.toLowerCase()
    if( !(genomeType in ['linear', 'circular']) ) {
        throw new IllegalArgumentException("Parameter --target_genome_type must be 'linear' or 'circular'. Received: ${params.target_genome_type}")
    }
    params.target_genome_type = genomeType

    params.threads = (params.threads ?: 1) as Integer
    if( params.threads < 1 ) {
        throw new IllegalArgumentException("Parameter --threads must be >= 1. Received: ${params.threads}")
    }

    def BQ_THR    = (params.bq_thr ?: 13).toString().toInteger()
    def DEPTH_THR = (params.depth_thr ?: 5).toString().toInteger()
    def QUAL_THR  = (params.qual_thr ?: 30).toString().toInteger()
    if( BQ_THR < 0 ) {
        throw new IllegalArgumentException("Parameter --bq_thr must be >= 0. Received: ${params.bq_thr}")
    }
    if( DEPTH_THR < 0 ) {
        throw new IllegalArgumentException("Parameter --depth_thr must be >= 0. Received: ${params.depth_thr}")
    }
    if( QUAL_THR < 0 ) {
        throw new IllegalArgumentException("Parameter --qual_thr must be >= 0. Received: ${params.qual_thr}")
    }
    if( params.cleantmpdatasize_bytes < 0 ) {
        throw new IllegalArgumentException("Parameter --cleantmpdatasize must be >= 0. Received: ${params.cleantmpdatasize}")
    }

    def targetReference = new File(params.target_reference.toString())
    if( !targetReference.exists() || !targetReference.isFile() ) {
        throw new IllegalArgumentException("Target reference FASTA was not found: ${params.target_reference}")
    }

    def requiredScripts = [
        'python/link_fastq.py',
        'python/draw_qc_plot.py',
        'python/check_human_ref.py',
        'python/format_fasta.py',
        'python/em_real_location.py',
        'python/normalise_PE_fastq.py',
        'python/format_pileup.py',
        'python/dedouble_pileup.py',
        'python/dedouble_depth.py',
        'python/mask_depth.py',
        'python/get_coverage.py',
        'python/coverage_plot_interactive.py',
        'python/extract_top_references.py',
        'python/extract_vcf_info.py',
        'python/built_nt_freq.py',
        'python/merge_tsv.py',
        'python/export_command_md.py',
        'python/cleanup_tmp_data.py',
    ]
    def missingScripts = requiredScripts.findAll { !new File(projectDir.toString(), it).exists() }
    if( missingScripts ) {
        throw new IllegalArgumentException("Required helper scripts are missing: ${missingScripts.join(', ')}")
    }

    def requiredEnvs = [
        'environment/seqtk_env.yml',
        'environment/plot_env.yml',
        'environment/fastp_env.yml',
        'environment/cutadapt_env.yml',
        'environment/bowtie2_env.yml',
        'environment/samtools_env.yml',
        'environment/spades_env.yml',
        'environment/minimap2_env.yml',
    ]
    def missingEnvs = requiredEnvs.findAll { !new File(projectDir.toString(), it).exists() }
    if( missingEnvs ) {
        throw new IllegalArgumentException("Required conda environment files are missing: ${missingEnvs.join(', ')}")
    }

    if( params.remove_human ) {
        params.human_db = resolvePathParam(params.human_db, 'human_db', false)
        def suffixes = ['.1.bt2', '.2.bt2', '.3.bt2', '.4.bt2', '.rev.1.bt2', '.rev.2.bt2']
        def missing = suffixes.findAll { !new File("${params.human_db}${it}").exists() }
        if( missing ) {
            throw new IllegalArgumentException("Human Bowtie2 index prefix is incomplete. Missing: ${missing.collect { "${params.human_db}${it}" }.join(', ')}")
        }
    }
}


// ========================================================================================
// Workflow Definition
// ========================================================================================

workflow {
    main:
    validateParams()
    def primerPaths = resolvePrimerFastaParams(params.fwd_primer_seq, params.rev_primer_seq)
    def discoveredSamples = discoverInputPairs(params.inputcsv.toString())
    params.project_name = resolveProjectNameFromOutdir(params.outdir)

    log.info "Discovered ${discoveredSamples.size()} paired-end sample(s) from ${params.inputcsv}"
    discoveredSamples.each { sample ->
        log.info " - ${sample.id} (${sample.sample_name}): ${sample.r1_original_name}, ${sample.r2_original_name}"
    }

    // Step 1: build the primary input channel.
    // Output: tuples containing sample metadata plus the raw R1/R2 FASTQ paths.
    Channel
        .fromList(discoveredSamples)
        .map { sample ->
            tuple(
                sample.subMap('id', 'sample_name', 'r1_canonical_name', 'r2_canonical_name', 'r1_original_name', 'r2_original_name'),
                file(sample.r1_path),
                file(sample.r2_path)
            )
        }
        .set { input_samples_ch }

    Channel
        .value(tuple(file(primerPaths.fwd_primer_seq), file(primerPaths.rev_primer_seq)))
        .set { primer_files_ch }

    // Step 2: hand the normalized input channel to the production subworkflow.
    // Output: the full MLPS paired-end processing graph under the chosen outdir.
    MLPS_PE_PIPELINE(input_samples_ch, primer_files_ch)

    onComplete:
    mergePublishedSummaryResults()
    exportPipelineCommandAudit()
    cleanupTemporaryTaskData()
    printCompletionSummary()
}
