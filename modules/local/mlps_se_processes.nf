// Create canonical paired FASTQ filenames inside the pipeline and publish them for traceability.
// Output: `${sample}.R1.fastq.gz` and `${sample}.R2.fastq.gz` under `fastq_links/`.
process NORMALIZE_INPUT_FASTQ {
    tag { meta.id }
    label 'light'
    publishDir "${params.outdir}/00_fastq_links", mode: 'symlink'
    conda "${projectDir}/environment/seqtk_env.yml"



    input:
    tuple val(meta), path(raw_r1, stageAs: 'input.R1.fastq.gz'), path(raw_r2, stageAs: 'input.R2.fastq.gz')

    output:
    tuple val(meta), path("${meta.id}.R1.fastq.gz"), path("${meta.id}.R2.fastq.gz"), emit: pair

    script:
    """
    ln -s "${raw_r1}" "${meta.id}.R1.fastq.gz"
    ln -s "${raw_r2}" "${meta.id}.R2.fastq.gz"
    """

    stub:
    """
    ln -s "${raw_r1}" "${meta.id}.R1.fastq.gz"
    ln -s "${raw_r2}" "${meta.id}.R2.fastq.gz"
    """
}

// Count reads and bases in a FASTQ/FASTA-like input and emit a one-line TSV summary.
// Output: `${sample}.${stage_name}.tsv` under `summary/`.
process COUNT_READS {
    tag { "${meta.id}:${stage_name}" }
    label 'light'
    conda "${projectDir}/environment/seqtk_env.yml"



    input:
    tuple val(meta), path(reads), val(stage_name)

    output:
    tuple val(meta), val(stage_name), path('*.tsv'), emit: count_file

    script:
    """
    count=\$(seqtk size "${reads}" 2>/dev/null | cut -f1 || echo 0)
    count=\${count:-0}
    totalbase=\$(seqtk size "${reads}" 2>/dev/null | cut -f2 || echo 0)
    totalbase=\${totalbase:-0}
    printf "%s\\t%s\\t%s\\n" "${meta.id}" "\${count}" "\${totalbase}" > "${meta.id}.${stage_name}.tsv"
    """

    stub:
    """
    printf "%s\\t1\\t12\\n" "${meta.id}" > "${meta.id}.${stage_name}.tsv"
    """
}

// Count paired reads/bases and emit one merged sample-level TSV summary row.
// Output: `${sample}.${stage_name}.tsv` under `summary/`.
process COUNT_PAIRED_READS {
    tag { "${meta.id}:${stage_name}" }
    label 'light'
    conda "${projectDir}/environment/seqtk_env.yml"



    input:
    tuple val(meta), path(read1), path(read2), val(stage_name)

    output:
    tuple val(meta), val(stage_name), path('*.tsv'), emit: count_file

    script:
    """
    count_r1=\$(seqtk size "${read1}" 2>/dev/null | cut -f1 || echo 0)
    count_r1=\${count_r1:-0}
    bases_r1=\$(seqtk size "${read1}" 2>/dev/null | cut -f2 || echo 0)
    bases_r1=\${bases_r1:-0}
    count_r2=\$(seqtk size "${read2}" 2>/dev/null | cut -f1 || echo 0)
    count_r2=\${count_r2:-0}
    bases_r2=\$(seqtk size "${read2}" 2>/dev/null | cut -f2 || echo 0)
    bases_r2=\${bases_r2:-0}
    printf "%s\\t%s\\t%s\\n" "${meta.id}" "\$((count_r1 + count_r2))" "\$((bases_r1 + bases_r2))" > "${meta.id}.${stage_name}.tsv"
    """

    stub:
    """
    printf "%s\\t2\\t24\\n" "${meta.id}" > "${meta.id}.${stage_name}.tsv"
    """
}

// Generate seqtk fqchk tables and PDF QC plots for each mate in a sample.
// Output: `PreQC/*.preQC.tsv` and `PreQC_PDF/*.preQC.pdf`.
process GENERATE_PREQC {
    tag { meta.id }
    publishDir "${params.outdir}/01_fastq_reports", mode: 'copy'
    conda "${projectDir}/environment/plot_env.yml"



    input:
    tuple val(meta), path(read1), path(read2)

    output:
    tuple val(meta), path('PreQC/*.preQC.tsv'), emit: report
    tuple val(meta), path('PreQC_PDF/*.preQC.pdf'), emit: pdf

    script:
    """
    mkdir -p PreQC PreQC_PDF

    seqtk fqchk "${read1}" > "PreQC/${meta.id}.R1.preQC.tsv"
    seqtk fqchk "${read2}" > "PreQC/${meta.id}.R2.preQC.tsv"
    python "${projectDir}/python/draw_qc_plot.py" \
        --input "PreQC/${meta.id}.R1.preQC.tsv" \
        --output "PreQC_PDF/${meta.id}.R1.preQC.pdf"
    python "${projectDir}/python/draw_qc_plot.py" \
        --input "PreQC/${meta.id}.R2.preQC.tsv" \
        --output "PreQC_PDF/${meta.id}.R2.preQC.pdf"
    """

    stub:
    """
    mkdir -p PreQC PreQC_PDF
    cat <<'EOF' > "PreQC/${meta.id}.R1.preQC.tsv"
min_len: 10; max_len: 10; avg_len: 10; 41 distinct quality values
POS     #bases  %A  %C  %G  %T  %N  avgQ  errQ  %low  %high
ALL     10      25  25  25  25  0   35    35    0     100
1       1       25  25  25  25  0   35    35    0     100
EOF
    cp "PreQC/${meta.id}.R1.preQC.tsv" "PreQC/${meta.id}.R2.preQC.tsv"
    touch "PreQC_PDF/${meta.id}.R1.preQC.pdf"
    touch "PreQC_PDF/${meta.id}.R2.preQC.pdf"
    """
}

// Apply paired-end quality filtering with fastp.
// Output: `${sample}.R1.f1.fastq.gz` and `${sample}.R2.f1.fastq.gz` under `Quality_trimming/`.
process QUALITY_TRIM {
    tag { meta.id }
    publishDir "${params.outdir}/02_quality_trimming", mode: 'symlink'
    conda "${projectDir}/environment/fastp_env.yml"



    input:
    tuple val(meta), path(read1), path(read2)

    output:
    tuple val(meta), path("${meta.id}.R1.f1.fastq.gz"), path("${meta.id}.R2.f1.fastq.gz"), emit: pair

    script:
    """
    fastp \
        -i "${read1}" \
        -I "${read2}" \
        -o "${meta.id}.R1.f1.fastq.gz" \
        -O "${meta.id}.R2.f1.fastq.gz" \
        --detect_adapter_for_pe \
        --disable_adapter_trimming \
        --cut_front false \
        --cut_tail false \
        --cut_right false \
        --qualified_quality_phred 20 \
        --unqualified_percent_limit 40 \
        --length_required 30 \
        --thread ${task.cpus}
    """

    stub:
    """
    cp "${read1}" "${meta.id}.R1.f1.fastq.gz"
    cp "${read2}" "${meta.id}.R2.f1.fastq.gz"
    """
}

// Export primary-alignment CIGAR strings and summarize their operation totals.
// Output: `${sample}.primarycigar.txt` in the stage folder and a helper TSV row for summary merging.
process SUMMARISE_PRIMARY_CIGAR {
    tag { "${meta.id}:${stage_name}" }
    publishDir({ "${params.outdir}/${stage_name}" }, mode: 'copy', pattern: '*.primarycigar.txt')
    conda "${projectDir}/environment/samtools_env.yml"



    input:
    tuple val(meta), path(bam), val(stage_name), val(prefix)

    output:
    tuple val(meta), path('*.primarycigar.txt'), emit: cigar_text
    tuple val(meta), path('*.primarycigar_row.tsv'), emit: cigar_row

    script:
    """
    samtools view -F 1024 "${bam}" | cut -f6 > "${meta.id}.primarycigar.txt"
    python "${projectDir}/python/cigar_summary.py" \
        --input "${meta.id}.primarycigar.txt" \
        --outputtsv "${meta.id}.primarycigar.tsv" \
        --prefix "${prefix}" \
        --samplename "${meta.id}"
    tail -n +2 "${meta.id}.primarycigar.tsv" > "${meta.id}.primarycigar_row.tsv"
    """

    stub:
    """
    cat <<'EOF' > "${meta.id}.primarycigar.txt"
10M
2S8M
EOF
    cat <<EOF > "${meta.id}.primarycigar.tsv"
sample\t${prefix}_hardclip_bases\t${prefix}_softclip_bases\t${prefix}_match_bases\t${prefix}_insertion_bases\t${prefix}_deletion_bases
${meta.id}\t0\t2\t18\t0\t0
EOF
    printf "%s\\t0\\t2\\t18\\t0\\t0\\n" "${meta.id}" > "${meta.id}.primarycigar_row.tsv"
    """
}

// Downsample reads to a fixed depth for adapter trimming and pilot compatibility.
// Output: `${sample}.R1.f2.fastq.gz` and `${sample}.R2.f2.fastq.gz` under `Quality_trimming/`.
process SUBSAMPLE_READS {
    tag { meta.id }
    publishDir "${params.outdir}/02_quality_trimming", mode: 'symlink'
    conda "${projectDir}/environment/seqtk_env.yml"



    input:
    tuple val(meta), path(read1), path(read2)

    output:
    tuple val(meta), path("${meta.id}.R1.f2.fastq.gz"), path("${meta.id}.R2.f2.fastq.gz"), emit: pair

    script:
    """
    seqtk sample -s100 "${read1}" 1000000 | gzip -c > "${meta.id}.R1.f2.fastq.gz"
    seqtk sample -s100 "${read2}" 1000000 | gzip -c > "${meta.id}.R2.f2.fastq.gz"
    """

    stub:
    """
    cp "${read1}" "${meta.id}.R1.f2.fastq.gz"
    cp "${read2}" "${meta.id}.R2.f2.fastq.gz"
    """
}

// Remove known adapter sequences from the subsampled paired reads.
// Output: `${sample}.R1.f3.fastq.gz` and `${sample}.R2.f3.fastq.gz` under `Adapter_trimming/`.
process ADAPTER_TRIM {
    tag { meta.id }
    publishDir "${params.outdir}/03_adapter_trimming", mode: 'symlink'
    conda "${projectDir}/environment/cutadapt_env.yml"



    input:
    tuple val(meta), path(read1), path(read2)

    output:
    tuple val(meta), path("${meta.id}.R1.f3.fastq.gz"), path("${meta.id}.R2.f3.fastq.gz"), emit: pair

    script:
    """
    cutadapt \
        -g "^file:${projectDir}/reference/known_adapters.fasta" \
        -G "^file:${projectDir}/reference/known_adapters.fasta" \
        -o "${meta.id}.R1.f3.fastq.gz" \
        -p "${meta.id}.R2.f3.fastq.gz" \
        "${read1}" "${read2}" \
        -e 0.1 \
        -m 15 \
        --quiet
    """

    stub:
    """
    cp "${read1}" "${meta.id}.R1.f3.fastq.gz"
    cp "${read2}" "${meta.id}.R2.f3.fastq.gz"
    """
}

// Validate that the configured Bowtie2 human reference is usable before depletion starts.
// Output: sentinel file `human_reference.ok` under `reference/HumanGenome/`.
process CHECK_HUMAN_REFERENCE {
    label 'light'
    publishDir "${params.outdir}/04_reference/HumanGenome", mode: 'copy', overwrite: false
    conda "${projectDir}/environment/bowtie2_env.yml"



    input:
    val(human_db)

    output:
    path("human_reference.ok"), emit: check

    script:
    """
    python "${projectDir}/python/check_human_ref.py" \
        --humanref "${human_db}" \
        --rebuild no \
        --newrefdir "${params.outdir}/04_reference/HumanGenome"
    touch human_reference.ok
    """

    stub:
    """
    touch human_reference.ok
    """
}

// Filter out reads mapping to the human reference and keep only non-human read pairs.
// Output: `${sample}.R1.f4.fastq.gz` and `${sample}.R2.f4.fastq.gz`.
process REMOVE_HUMAN_READS {
    tag { meta.id }
    publishDir "${params.outdir}/05_nonhuman", mode: 'symlink'
    conda "${projectDir}/environment/bowtie2_env.yml"



    input:
    tuple val(meta), path(read1), path(read2)

    output:
    tuple val(meta), path("${meta.id}.R1.f4.fastq.gz"), path("${meta.id}.R2.f4.fastq.gz"), emit: pair

    script:
    """
    bowtie2 \
        -x "${params.human_db}" \
        -q \
        --time \
        --reorder \
        -p ${task.cpus} \
        -1 "${read1}" \
        -2 "${read2}" \
        2> "${meta.id}.bowtie2.log" \
        | samtools view -h -f13 -F3840 - \
        | samtools sort -n - \
        | samtools fastq -n \
            -1 "${meta.id}.R1.f4.fastq.gz" \
            -2 "${meta.id}.R2.f4.fastq.gz" \
            -0 /dev/null \
            -s /dev/null
    """

    stub:
    """
    cp "${read1}" "${meta.id}.R1.f4.fastq.gz"
    cp "${read2}" "${meta.id}.R2.f4.fastq.gz"
    """
}

// Build cleaned linear and doubled target references used by downstream alignment steps.
// Output: `ref.format.fasta` and `ref.double.fasta` under `reference/`.
process FORMAT_TARGET_REFERENCE {
    label 'light'
    publishDir "${params.outdir}/04_reference", mode: 'copy'
    conda "${projectDir}/environment/samtools_env.yml"



    input:
    path(target_reference)

    output:
    tuple path("ref.format.fasta"), path("ref.double.fasta"), emit: refs
    path("ref.format.fasta"), emit: formatted
    path("ref.double.fasta"), emit: doubled

    script:
    """
    python "${projectDir}/python/format_fasta.py" "${target_reference}" \
        --uppercase yes \
        --removegap yes \
        --double no \
        --output ref.format.fasta

    python "${projectDir}/python/format_fasta.py" "${target_reference}" \
        --uppercase yes \
        --removegap yes \
        --double yes \
        --output ref.double.fasta
    """

    stub:
    """
    cp "${target_reference}" ref.format.fasta
    cp "${target_reference}" ref.double.fasta
    """
}

// Align non-human read pairs to the formatted reference set.
// Output: `${sample}.mixgenome.bam` plus its index in the work directory.
process MIXGENOME_ALIGN {
    tag { meta.id }
    publishDir "${params.outdir}/06_mixgenome", mode: 'symlink'
    conda "${projectDir}/environment/samtools_env.yml"



    input:
    tuple val(meta), path(read1), path(read2), path(ref_format), path(ref_double)

    output:
    tuple val(meta), path('*.mixgenome.bam'), emit: bam

    script:
    def mappingRef = params.target_genome_type == 'circular' ? ref_double.getName() : ref_format.getName()
    """
    bwa index "${mappingRef}" >/dev/null 2>&1
    bwa mem -a "${mappingRef}" "${read1}" "${read2}" \
        | samtools view -h -F3588 - \
        | samtools sort -@ ${task.cpus} -o "${meta.id}.mixgenome.bam" -
    samtools index "${meta.id}.mixgenome.bam"
    """

    stub:
    """
    touch "${meta.id}.mixgenome.bam"
    """
}

// Reallocate multi-mapped reads across references with the EM helper script.
// Output: per-read assignments TSV and per-reference abundance summary TSV.
process EM_REALLOCATE_REFERENCES {
    tag { meta.id }
    publishDir "${params.outdir}/06_mixgenome", mode: 'symlink', pattern: '*.em_assignments.tsv'
    publishDir "${params.outdir}/06_mixgenome", mode: 'copy', pattern: '*.em_assignments_ref_summary.tsv'
    conda "${projectDir}/environment/samtools_env.yml"



    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path('*.em_assignments.tsv'), emit: assignments
    tuple val(meta), path('*.em_assignments_ref_summary.tsv'), emit: summary

    script:
    """
    python "${projectDir}/python/em_real_location.py" \
        --bam "${bam}" \
        --out "${meta.id}.em_assignments.tsv"
    """

    stub:
    """
    cat <<'EOF' > "${meta.id}.em_assignments.tsv"
read_id\tref\tposterior
r1\tref1\t1.0
EOF
    cat <<'EOF' > "${meta.id}.em_assignments_ref_summary.tsv"
Reference\tTotal_posterior_mass\tRelative_abundance
ref1\t1.0\t1.0
EOF
    """
}

// Select the top-supported reference FASTA for each sample.
// Output: `${sample}.selected.fasta` under `reference/`.
process EXTRACT_TOP_REFERENCE {
    tag { meta.id }
    publishDir "${params.outdir}/04_reference", mode: 'copy'
    conda "${projectDir}/environment/samtools_env.yml"



    input:
    tuple val(meta), path(ref_summary), path(formatted_reference)

    output:
    tuple val(meta), path('*.selected.fasta'), emit: reference

    script:
    """
    python "${projectDir}/python/extract_top_references.py" \
        --tsv "${ref_summary}" \
        --fasta "${formatted_reference}" \
        --out "${meta.id}.selected.fasta" \
        --topn 1
    """

    stub:
    """
    cat <<'EOF' > "${meta.id}.selected.fasta"
>${meta.id}
ACGTACGTACGT
EOF
    """
}

// Prepare per-sample linear and doubled reference FASTA files for consensus mapping.
// Output: `${sample}.selected.fasta` and `${sample}.selected.double.fasta`.
process PREPARE_SAMPLE_REFERENCE {
    tag { meta.id }
    publishDir "${params.outdir}/04_reference", mode: 'copy'
    conda "${projectDir}/environment/samtools_env.yml"



    input:
    tuple val(meta), path(reference_fasta, stageAs: 'input_reference.fasta')

    output:
    tuple val(meta), path('*.selected.fasta'), path('*.selected.double.fasta'), emit: refs

    script:
    """
    cp "${reference_fasta}" "${meta.id}.selected.fasta"
    python "${projectDir}/python/format_fasta.py" "${reference_fasta}" \
        --uppercase yes \
        --removegap yes \
        --double yes \
        --output "${meta.id}.selected.double.fasta"
    """

    stub:
    """
    cp "${reference_fasta}" "${meta.id}.selected.fasta"
    cp "${reference_fasta}" "${meta.id}.selected.double.fasta"
    """
}

// Extract target-aligned read pairs from the mixed-reference BAM.
// Output: `${sample}.R1.target.fastq.gz` and `${sample}.R2.target.fastq.gz` under `mixgenome/`.
process BUILD_TARGET_READ_POOL {
    tag { meta.id }
    publishDir "${params.outdir}/06_mixgenome", mode: 'symlink'
    conda "${projectDir}/environment/samtools_env.yml"



    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${meta.id}.R1.target.fastq.gz"), path("${meta.id}.R2.target.fastq.gz"), emit: pair

    script:
    """
    samtools fastq -f1 -F256 -n "${bam}" \
        -1 "${meta.id}.R1.target.fastq.gz" \
        -2 "${meta.id}.R2.target.fastq.gz" \
        -0 /dev/null \
        -s /dev/null
    """

    stub:
    """
    cat <<'EOF' | gzip -c > "${meta.id}.R1.target.fastq.gz"
@${meta.id}_1
ACGTACGTACGT
+
FFFFFFFFFFFF
EOF
    cat <<'EOF' | gzip -c > "${meta.id}.R2.target.fastq.gz"
@${meta.id}_1
TGCATGCATGCA
+
FFFFFFFFFFFF
EOF
    """
}

// Build the first-pass alignment against the chosen sample reference.
// Output: `${sample}.L1.bam` plus its index in the work directory.
process L1_ALIGN {
    tag { meta.id }
    publishDir "${params.outdir}/07_L1_results", mode: 'copy'
    conda "${projectDir}/environment/samtools_env.yml"



    input:
    tuple val(meta), path(read1), path(read2), path(reference_fasta), path(reference_double)

    output:
    tuple val(meta), path('*.L1.bam'), emit: bam

    script:
    def mappingRef = params.target_genome_type == 'circular' ? reference_double.getName() : reference_fasta.getName()
    def dedupCmd = params.dedup
        ? """bwa mem "${mappingRef}" "${read1}" "${read2}" | samtools view -h -f1 -F3844 - | samtools sort -@ ${task.cpus} -n | samtools fixmate -m - - | samtools sort -@ ${task.cpus} | samtools markdup - "${meta.id}.L1.bam" """
        : """bwa mem "${mappingRef}" "${read1}" "${read2}" | samtools view -h -f1 -F3844 - | samtools sort -@ ${task.cpus} -o "${meta.id}.L1.bam" -"""
    """
    bwa index "${mappingRef}" >/dev/null 2>&1
    ${dedupCmd}
    samtools index "${meta.id}.L1.bam"
    """

    stub:
    """
    touch "${meta.id}.L1.bam"
    """
}

// Summarize key SAM flag counts for each L1 alignment BAM.
// Output: `${sample}.flags.tsv` under `summary/`.
process COUNT_ALIGNMENT_FLAGS {
    tag { meta.id }
    label 'light'
    conda "${projectDir}/environment/samtools_env.yml"



    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path('*.flags.tsv'), emit: flag_file

    script:
    """
    samtools view "${bam}" | \
    awk -v sample="${meta.id}" '
    {
        flag=\$2
        rev=(int(flag/16)%2)
        dup=(int(flag/1024)%2)
        if( rev && dup ) rev_dup++
        else if( rev ) rev_only++
        else if( dup ) fwd_dup++
        else fwd_only++
    }
    END{
        printf "%s\\t%d\\t%d\\t%d\\t%d\\n", sample, fwd_only+0, rev_only+0, fwd_dup+0, rev_dup+0
    }' > "${meta.id}.flags.tsv"
    """

    stub:
    """
    printf "%s\\t0\\t0\\t0\\t0\\n" "${meta.id}" > "${meta.id}.flags.tsv"
    """
}

// Convert an alignment BAM into a pileup and preserve the original reference alongside it.
// Output: `${sample}.pileup` and `${sample}.original_reference.fasta`.
process BUILD_PILEUP {
    tag { meta.id }
    publishDir({ bam.getName().endsWith('.L1.bam') ? "${params.outdir}/07_L1_results" : bam.getName().endsWith('.L2.bam') ? "${params.outdir}/09_L2_results" : "${params.outdir}/08_assembly" }, mode: 'copy')
    conda "${projectDir}/environment/samtools_env.yml"



    input:
    tuple val(meta), path(bam), path(reference_fasta), path(reference_double)

    output:
    tuple val(meta), path('*.pileup'), path('*.original_reference.fasta'), emit: pileup

    script:
    def mappingRef = params.target_genome_type == 'circular' ? reference_double.getName() : reference_fasta.getName()
    def pileupSource = params.target_genome_type == 'circular' ? "${meta.id}.pileup.double" : "${meta.id}.pileup"
    def dedoubleCmd = params.target_genome_type == 'circular'
        ? """python "${projectDir}/python/dedouble_pileup.py" --input "${meta.id}.pileup.double" --originalref "${reference_fasta.getName()}" --output "${meta.id}.pileup" """
        : 'true'
    """
    cp "${reference_fasta}" "${meta.id}.original_reference.fasta"
    samtools faidx "${mappingRef}"
    samtools mpileup -aa -A -B -Q 0 -q 0 -f "${mappingRef}" "${bam}" > "${pileupSource}"
    ${dedoubleCmd}
    """

    stub:
    """
    cp "${reference_fasta}" "${meta.id}.original_reference.fasta"
    cat <<'EOF' > "${meta.id}.pileup"
ref\t1\tA\t1\t.\tF
ref\t2\tC\t1\t.\tF
EOF
    """
}

// Format a pileup into consensus, patched-reference, nt-frequency, and FASTA summary deliverables.
// Output: `*.ntfreq.csv`, `*.consensus.fasta`, `*.bases.fasta`, `*.patchedref.fasta`, `*.variantsignals.fasta`.
process FORMAT_PILEUP {
    tag { "${meta.id}:${stage_name}" }
    publishDir({ "${params.outdir}/${stage_name}" }, mode: 'copy')
    conda "${projectDir}/environment/samtools_env.yml"



    input:
    tuple val(meta), path(pileup, stageAs: 'input.pileup'), path(original_reference), val(stage_name)

    output:
    tuple val(meta), path('*.ntfreq.csv'), emit: ntfreq
    tuple val(meta), path('*.consensus.fasta'), emit: consensus_fasta
    tuple val(meta), path('*.bases.fasta'), emit: bases_fasta
    tuple val(meta), path('*.patchedref.fasta'), emit: patched_reference
    tuple val(meta), path('*.variantsignals.fasta'), emit: variant_signal
    tuple val(meta), path('*.fasta_summary_row.tsv'), emit: fasta_summary_row

    script:
    def summaryPrefix = stage_name.startsWith('07_') ? 'L1' : stage_name.startsWith('09_') ? 'L2' : 'assembly'
    def summaryName = stage_name.startsWith('07_')
        ? "${meta.id}.07_samtoolsseq_summmary.tsv"
        : stage_name.startsWith('09_')
            ? "${meta.id}.09_samtoolsseq_summmary.tsv"
            : null
    def summaryTarget = stage_name.startsWith('09_') ? "${meta.id}.consensus.fasta" : "${meta.id}.bases.fasta"
    def summaryCmd = summaryName
        ? """python "${projectDir}/python/get_fasta_summary.py" \
        --input "${summaryTarget}" \
        --prefix "${summaryPrefix}" \
        --outputtsv "${summaryName}" \
        --samplename "${meta.id}" \
        --seqtype NT
    tail -n +2 "${summaryName}" > "${meta.id}.fasta_summary_row.tsv" """
        : """printf "%s\\t0\\t0\\t0\\t0\\n" "${meta.id}" > "${meta.id}.fasta_summary_row.tsv" """
    """
    cp "${pileup}" "${meta.id}.pileup"
    python "${projectDir}/python/format_pileup.py" "${meta.id}.pileup"
    ${summaryCmd}
    """

    stub:
    """
    cat <<'EOF' > "${meta.id}.ntfreq.csv"
Ref,Pos,A,C,G,T,N,Depth,CasedConsensus,CalledConsensus,PatchedWithRef,VariantSignal
ref,1,1,0,0,0,0,1,A,A,A,A
EOF
    cat <<'EOF' > "${meta.id}.consensus.fasta"
>${meta.id}_caseconsensus
A
EOF
    cat <<'EOF' > "${meta.id}.bases.fasta"
>${meta.id}_calledconsensus
A
EOF
    cat <<'EOF' > "${meta.id}.patchedref.fasta"
>${meta.id}_patchedref
A
EOF
    cat <<'EOF' > "${meta.id}.variantsignals.fasta"
>${meta.id}_variantsignals
A
EOF
    printf "%s\\t1\\t0\\t0\\t0\\n" "${meta.id}" > "${meta.id}.fasta_summary_row.tsv"
    """
}

// Perform paired-end digital normalization before de novo assembly.
// Output: `${sample}.R1.normalised.fastq.gz` and `${sample}.R2.normalised.fastq.gz` under `assembly/`.
process NORMALISE_FOR_ASSEMBLY {
    tag { meta.id }
    publishDir "${params.outdir}/08_assembly", mode: 'symlink'
    conda "${projectDir}/environment/samtools_env.yml"



    input:
    tuple val(meta), path(read1), path(read2)

    output:
    tuple val(meta), path("${meta.id}.R1.normalised.fastq.gz"), path("${meta.id}.R2.normalised.fastq.gz"), emit: pair

    script:
    """
    python "${projectDir}/python/normalise_PE_fastq.py" \
        -1 "${read1}" \
        -2 "${read2}" \
        -o "${meta.id}.R1.normalised.fastq.gz" \
        -O "${meta.id}.R2.normalised.fastq.gz"
    """

    stub:
    """
    cp "${read1}" "${meta.id}.R1.normalised.fastq.gz"
    cp "${read2}" "${meta.id}.R2.normalised.fastq.gz"
    """
}

// Assemble normalized read pairs into contigs with SPAdes.
// Output: `${sample}.contigs.fasta` under `assembly/`.
process RUN_SPADES {
    tag { meta.id }
    publishDir "${params.outdir}/08_assembly", mode: 'copy'
    conda "${projectDir}/environment/spades_env.yml"



    input:
    tuple val(meta), path(read1), path(read2)

    output:
    tuple val(meta), path('*.contigs.fasta'), emit: contigs

    script:
    """
    spades.py \
        -1 "${read1}" \
        -2 "${read2}" \
        --only-assembler \
        --cov-cutoff 5 \
        -k 21 \
        -o "${meta.id}_spades_output"

    cp "${meta.id}_spades_output/contigs.fasta" "${meta.id}.contigs.fasta"
    """

    stub:
    """
    cat <<'EOF' > "${meta.id}.contigs.fasta"
>${meta.id}_contig1
ACGTACGTACGT
EOF
    """
}

// Align assembled contigs back to the sample reference.
// Output: `${sample}.assembly.bam` plus its index in the work directory.
process ALIGN_ASSEMBLY_TO_REFERENCE {
    tag { meta.id }
    publishDir "${params.outdir}/08_assembly", mode: 'symlink'
    conda "${projectDir}/environment/minimap2_env.yml"



    input:
    tuple val(meta), path(contigs), path(reference_fasta), path(reference_double)

    output:
    tuple val(meta), path('*.assembly.bam'), emit: bam

    script:
    def mappingRef = params.target_genome_type == 'circular' ? reference_double.getName() : reference_fasta.getName()
    """
    minimap2 -ax map-ont "${mappingRef}" "${contigs}" \
        | samtools view -h -F3844 - \
        | samtools sort -@ ${task.cpus} -o "${meta.id}.assembly.bam" -
    samtools index "${meta.id}.assembly.bam"
    """

    stub:
    """
    touch "${meta.id}.assembly.bam"
    """
}

// Prepare the patched reference pair used for the L2 remapping stage.
// Output: `${sample}.patchedref.original.fasta` and `${sample}.patchedref.fasta`.
process PREPARE_PATCHED_REFERENCE {
    tag { meta.id }
    publishDir "${params.outdir}/08_assembly", mode: 'copy', pattern: '*.patchedref.original.fasta'
    publishDir "${params.outdir}/04_reference", mode: 'copy', pattern: '*.patchedref.fasta'
    conda "${projectDir}/environment/samtools_env.yml"



    input:
    tuple val(meta), path(patched_reference, stageAs: 'input_patchedref.fasta')

    output:
    tuple val(meta), path('*.patchedref.original.fasta'), path('*.patchedref.fasta'), emit: refs

    script:
    def prepCmd = params.target_genome_type == 'circular'
        ? """python "${projectDir}/python/format_fasta.py" "${patched_reference}" --uppercase yes --removegap yes --double yes --output "${meta.id}.patchedref.fasta" """
        : """cp "${patched_reference}" "${meta.id}.patchedref.fasta" """
    """
    cp "${patched_reference}" "${meta.id}.patchedref.original.fasta"
    ${prepCmd}
    """

    stub:
    """
    cp "${patched_reference}" "${meta.id}.patchedref.original.fasta"
    cp "${patched_reference}" "${meta.id}.patchedref.fasta"
    """
}

// Remap target reads to the patched reference to build the final consensus stage.
// Output: `${sample}.L2.bam` plus its index in the work directory.
process L2_ALIGN {
    tag { meta.id }
    publishDir "${params.outdir}/09_L2_results", mode: 'copy'
    conda "${projectDir}/environment/samtools_env.yml"



    input:
    tuple val(meta), path(read1), path(read2), path(original_reference), path(mapping_reference)

    output:
    tuple val(meta), path('*.L2.bam'), emit: bam

    script:
    def dedupCmd = params.dedup
        ? """bwa mem "${mapping_reference.getName()}" "${read1}" "${read2}" | samtools view -h -f1 -F3844 - | samtools sort -@ ${task.cpus} -n | samtools fixmate -m - - | samtools sort -@ ${task.cpus} | samtools markdup - "${meta.id}.L2.bam" """
        : """bwa mem "${mapping_reference.getName()}" "${read1}" "${read2}" | samtools view -h -f1 -F3844 - | samtools sort -@ ${task.cpus} -o "${meta.id}.L2.bam" -"""
    """
    bwa index "${mapping_reference.getName()}" >/dev/null 2>&1
    ${dedupCmd}
    samtools index "${meta.id}.L2.bam"
    """

    stub:
    """
    touch "${meta.id}.L2.bam"
    """
}

// Expand nucleotide frequency tables into the downstream non-aligned nucleotide report format.
// Output: `${sample}.L2nt.csv` under `10_nonaligned_nt_counts/`.
process BUILD_NT_COUNTS {
    tag { meta.id }
    publishDir "${params.outdir}/10_nonaligned_nt_counts", mode: 'copy'
    conda "${projectDir}/environment/samtools_env.yml"



    input:
    tuple val(meta), path(ntfreq_csv)

    output:
    tuple val(meta), path('*.L2nt.csv'), emit: nt_csv

    script:
    """
    python "${projectDir}/python/built_nt_freq.py" "${ntfreq_csv}" "${meta.id}.L2nt.csv"
    """

    stub:
    """
    cat <<'EOF' > "${meta.id}.L2nt.csv"
non-aligned_pos,consensus_nt,nt,nt_count,nt_percent
1,A,A,1,100.00
EOF
    """
}

process COUNT_Q40_PRIMARY_READS {
    tag { meta.id }
    label 'light'
    conda "${projectDir}/environment/samtools_env.yml"



    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path('*.multi_q40.tsv'), emit: count_file

    script:
    """
    count=\$(samtools view -F 256 -q 40 "${bam}" -c 2>/dev/null || echo 0)
    printf "%s\\t%s\\n" "${meta.id}" "\${count}" > "${meta.id}.multi_q40.tsv"
    """

    stub:
    """
    printf "%s\\t1\\n" "${meta.id}" > "${meta.id}.multi_q40.tsv"
    """
}

process WRITE_BEST_REFERENCE {
    tag { meta.id }
    label 'light'
    conda "${projectDir}/environment/samtools_env.yml"



    input:
    tuple val(meta), path(selected_reference)

    output:
    tuple val(meta), path('*.best_ref.tsv'), emit: summary_file

    script:
    """
    seqname=\$(grep -m1 '^>' "${selected_reference}" 2>/dev/null | cut -d' ' -f1 | sed 's/^>//' || echo 'Failed')
    printf "%s\\t%s\\n" "${meta.id}" "\${seqname:-Failed}" > "${meta.id}.best_ref.tsv"
    """

    stub:
    """
    printf "%s\\t%s\\n" "${meta.id}" "${meta.id}" > "${meta.id}.best_ref.tsv"
    """
}

// Build coverage depth plus a mask BED.
// Linear mode uses the raw samtools depth output directly, while circular mode
// first dedoubles the depth table against the stage's original reference fasta.
process DEPTH_AND_MASK {
    tag { "${meta.id}:${stage_name}" }
    publishDir({ "${params.outdir}/${stage_name}" }, mode: 'copy')
    conda "${projectDir}/environment/samtools_env.yml"



    input:
    tuple val(meta), path(bam), path(original_reference), val(stage_name)

    output:
    tuple val(meta), path('*.depth.txt'), path('*.L*mask.bed'), emit: depth_mask

    script:
    def maskName = stage_name.startsWith('07_') ? "${meta.id}.L1mask.bed" : stage_name.startsWith('09_') ? "${meta.id}.L2mask.bed" : "${meta.id}.mask.bed"
    def depthCmd = params.target_genome_type == 'circular'
        ? """samtools depth -a "${bam}" > "${meta.id}.depth.double"
    python "${projectDir}/python/dedouble_depth.py" \
        --input "${meta.id}.depth.double" \
        --originalref "${original_reference}" \
        --output "${meta.id}.depth.txt" """
        : """samtools depth -a "${bam}" > "${meta.id}.depth.txt" """
    """
    ${depthCmd}
    python "${projectDir}/python/mask_depth.py" \
        --input "${meta.id}.depth.txt" \
        --bq ${params.depth_thr} \
        --output "${maskName}"
    """

    stub:
    def maskName = stage_name.startsWith('07_') ? "${meta.id}.L1mask.bed" : stage_name.startsWith('09_') ? "${meta.id}.L2mask.bed" : "${meta.id}.mask.bed"
    """
    cat <<'EOF' > "${meta.id}.depth.txt"
ref\t1\t10
ref\t2\t12
EOF
    : > "${maskName}"
    """
}

process SUMMARISE_COVERAGE {
    tag { "${meta.id}:${stage_name}" }
    label 'light'
    conda "${projectDir}/environment/plot_env.yml"



    input:
    tuple val(meta), path(depth_file), val(stage_name)

    output:
    tuple val(meta), path('*.coverage_row.tsv'), emit: coverage_row

    script:
    """
    python "${projectDir}/python/get_coverage.py" \
        --input "${depth_file}" \
        --prefix "${meta.id}" \
        --output "${meta.id}.${stage_name}.coverage.tsv"
    tail -n +2 "${meta.id}.${stage_name}.coverage.tsv" > "${meta.id}.${stage_name}.coverage_row.tsv"
    """

    stub:
    """
    printf "%s\\tref:10.000000\\tref:10.000000\\t10.000000\\t10.000000\\n" "${meta.id}" > "${meta.id}.${stage_name}.coverage_row.tsv"
    """
}

process BUILD_COVERAGE_PLOT {
    label 'light'
    publishDir "${params.outdir}/summary", mode: 'copy'
    conda "${projectDir}/environment/plot_env.yml"



    input:
    path(depth_files)
    val(output_name)

    output:
    path(output_name), emit: html

    script:
    def files = depth_files.collect { "\"${it.getName()}\"" }.join(' ')
    """
    mkdir -p depth_inputs
    cp ${files} depth_inputs/
    python "${projectDir}/python/coverage_plot_interactive.py" \
        --indir depth_inputs \
        --yaxis log10 \
        --output "${output_name}"
    """

    stub:
    """
    touch "${output_name}"
    """
}

process BCFTOOLS_VARIANTS {
    tag { "${meta.id}:${stage_name}" }
    publishDir(
        { "${params.outdir}/${stage_name}" },
        mode: 'copy',
        saveAs: { filename -> filename.endsWith('.variant_status.tsv') ? null : filename }
    )
    conda "${projectDir}/environment/samtools_env.yml"



    input:
    tuple val(meta), path(bam), path(reference_fasta), path(mask_bed), val(stage_name)

    output:
    tuple val(meta), path('*.bcfconsensus.fasta'), emit: consensus_fasta
    tuple val(meta), path('*.vcf.gz'), emit: raw_vcf
    tuple val(meta), path('*.norm.vcf.gz'), emit: norm_vcf
    tuple val(meta), path('*.norm.vcf.tsv'), emit: norm_vcf_tsv
    tuple val(meta), path('*.filter.vcf.gz'), emit: filter_vcf
    tuple val(meta), path('*.variant_status.tsv'), emit: status_file
    tuple val(meta), path('*.bcfseq_summary_row.tsv'), emit: fasta_summary_row

    script:
    def summaryName = stage_name.startsWith('07_') ? "${meta.id}.07_bcfseq_summmary.tsv" : stage_name.startsWith('09_') ? "${meta.id}.09_bcfseq_summmary.tsv" : "${meta.id}.${stage_name}_bcfseq_summmary.tsv"
    def summaryPrefix = stage_name.startsWith('07_') ? 'L1' : stage_name.startsWith('09_') ? 'L2' : stage_name
    """
    bcftools mpileup -Ou -Q ${params.bq_thr} -f "${reference_fasta}" "${bam}" \
      | bcftools call -mv --ploidy 1 -Oz -o "${meta.id}.vcf.gz"
    bcftools index "${meta.id}.vcf.gz"
    bcftools norm -f "${reference_fasta}" -Oz -o "${meta.id}.norm.vcf.gz" "${meta.id}.vcf.gz"
    bcftools index "${meta.id}.norm.vcf.gz"
    python "${projectDir}/python/extract_vcf_info.py" \
      --input "${meta.id}.norm.vcf.gz" \
      --outputtsv "${meta.id}.norm.vcf.tsv"
    bcftools view -i "DP>${params.depth_thr} && QUAL>${params.qual_thr}" "${meta.id}.norm.vcf.gz" -Oz -o "${meta.id}.filter.vcf.gz"
    bcftools index "${meta.id}.filter.vcf.gz"
    bcftools consensus -f "${reference_fasta}" -m "${mask_bed}" "${meta.id}.filter.vcf.gz" \
      | sed "1s/^>.*/>${meta.id}_bcfconsensus/" > "${meta.id}.bcfconsensus.fasta"
    python "${projectDir}/python/get_fasta_summary.py" \
      --input "${meta.id}.bcfconsensus.fasta" \
      --prefix "${summaryPrefix}" \
      --outputtsv "${summaryName}" \
      --samplename "${meta.id}" \
      --seqtype NT
    tail -n +2 "${summaryName}" > "${meta.id}.bcfseq_summary_row.tsv"
    variant_sites=\$(bcftools view -H "${meta.id}.vcf.gz" 2>/dev/null | wc -l || echo 0)
    high_qual_variants=\$(bcftools view -H "${meta.id}.filter.vcf.gz" 2>/dev/null | wc -l || echo 0)
    printf "%s\\t%s\\t%s\\n" "${meta.id}" "\${variant_sites:-0}" "\${high_qual_variants:-0}" > "${meta.id}.variant_status.tsv"
    """

    stub:
    """
    touch "${meta.id}.vcf.gz" "${meta.id}.norm.vcf.gz" "${meta.id}.filter.vcf.gz"
    cat <<'EOF' > "${meta.id}.norm.vcf.tsv"
CHROM\tPOS\tREF\tALT\tQUAL\tDP\tQ13_ref_fwd\tQ13_ref_rev\tQ13_alt_fwd\tQ13_alt_rev\tallele_freq\terror_rate
ref\t1\tA\tG\t40\t10\t5\t0\t5\t0\t0.500000\t0.0001
EOF
    cat <<'EOF' > "${meta.id}.bcfconsensus.fasta"
>${meta.id}_bcfconsensus
A
EOF
    printf "%s\\t1\\t0\\t0\\t0\\n" "${meta.id}" > "${meta.id}.bcfseq_summary_row.tsv"
    printf "%s\\t1\\t1\\n" "${meta.id}" > "${meta.id}.variant_status.tsv"
    """
}

process WRITE_VARIANT_STATUS_PLACEHOLDER {
    tag { "${meta.id}:${stage_name}" }
    label 'light'
    conda "${projectDir}/environment/seqtk_env.yml"



    input:
    tuple val(meta), val(stage_name)

    output:
    tuple val(meta), path('*.variant_status.tsv'), emit: status_file

    script:
    """
    printf "%s\\tCircular-NA\\tCircular-NA\\n" "${meta.id}" > "${meta.id}.${stage_name}.variant_status.tsv"
    """

    stub:
    """
    printf "%s\\tCircular-NA\\tCircular-NA\\n" "${meta.id}" > "${meta.id}.${stage_name}.variant_status.tsv"
    """
}

// Merge per-sample TSV outputs into a single summary table with a provided header.
// Output: one summary TSV under `summary/`.
process MERGE_SUMMARY_TABLE {
    label 'light'
    publishDir "${params.outdir}/summary", mode: 'copy'
    conda "${projectDir}/environment/seqtk_env.yml"



    input:
    path(count_files)
    val(output_name)
    val(header_line)

    output:
    path(output_name), emit: summary

    script:
    def files = count_files.collect { "\"${it.getName()}\"" }.join(' ')
    """
    {
        printf "%s\\n" "${header_line}"
        cat ${files}
    } > "${output_name}"
    """

    stub:
    """
    printf "%s\\n" "${header_line}" > "${output_name}"
    """
}

// Concatenate FASTA collections into a combined summary FASTA file.
// Output: one merged FASTA under `summary/`.
process MERGE_FASTA_COLLECTION {
    label 'light'
    publishDir "${params.outdir}/summary", mode: 'copy'
    conda "${projectDir}/environment/seqtk_env.yml"



    input:
    path(fasta_files)
    val(output_name)

    output:
    path(output_name), emit: fasta

    script:
    def files = fasta_files.collect { "\"${it.getName()}\"" }.join(' ')
    """
    cat ${files} > "${output_name}"
    """

    stub:
    """
    touch "${output_name}"
    """
}
