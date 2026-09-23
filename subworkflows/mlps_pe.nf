include {
    NORMALIZE_INPUT_FASTQ;
    COUNT_PAIRED_READS as COUNT_PREQC_READS;
    COUNT_PAIRED_READS as COUNT_POSTQF_READS;
    COUNT_PAIRED_READS as COUNT_POSTSUB_READS;
    COUNT_PAIRED_READS as COUNT_POSTAT_READS;
    COUNT_PAIRED_READS as COUNT_NONHUMAN_READS;
    COUNT_PAIRED_READS as COUNT_TARGET_READS;
    COUNT_PAIRED_READS as COUNT_ASSEMBLY_READS;
    COUNT_READS as COUNT_CONTIG_READS;
    GENERATE_PREQC;
    QUALITY_TRIM;
    SUBSAMPLE_READS;
    ADAPTER_TRIM;
    CHECK_HUMAN_REFERENCE;
    REMOVE_HUMAN_READS;
    FORMAT_TARGET_REFERENCE;
    MIXGENOME_ALIGN;
    EM_REALLOCATE_REFERENCES;
    EXTRACT_TOP_REFERENCE;
    PREPARE_SAMPLE_REFERENCE;
    BUILD_TARGET_READ_POOL;
    COUNT_Q40_PRIMARY_READS;
    WRITE_BEST_REFERENCE;
    L1_ALIGN;
    COUNT_ALIGNMENT_FLAGS as COUNT_L1_ALIGNMENT_FLAGS;
    COUNT_ALIGNMENT_FLAGS as COUNT_L2_ALIGNMENT_FLAGS;
    SUMMARISE_PRIMARY_CIGAR as SUMMARISE_L1_PRIMARY_CIGAR;
    SUMMARISE_PRIMARY_CIGAR as SUMMARISE_L2_PRIMARY_CIGAR;
    DEPTH_AND_MASK as BUILD_L1_DEPTH_MASK;
    DEPTH_AND_MASK as BUILD_L2_DEPTH_MASK;
    SUMMARISE_COVERAGE as SUMMARISE_L1_COVERAGE;
    SUMMARISE_COVERAGE as SUMMARISE_L2_COVERAGE;
    BUILD_COVERAGE_PLOT as BUILD_L1_COVERAGE_PLOT;
    BUILD_COVERAGE_PLOT as BUILD_L2_COVERAGE_PLOT;
    BCFTOOLS_VARIANTS as BUILD_L1_VARIANTS;
    BCFTOOLS_VARIANTS as BUILD_L2_VARIANTS;
    WRITE_VARIANT_STATUS_PLACEHOLDER;
    BUILD_PILEUP as BUILD_L1_PILEUP;
    BUILD_PILEUP as BUILD_ASSEMBLY_PILEUP;
    BUILD_PILEUP as BUILD_L2_PILEUP;
    FORMAT_PILEUP as FORMAT_L1_PILEUP;
    FORMAT_PILEUP as FORMAT_ASSEMBLY_PILEUP;
    FORMAT_PILEUP as FORMAT_L2_PILEUP;
    NORMALISE_FOR_ASSEMBLY;
    RUN_SPADES;
    ALIGN_ASSEMBLY_TO_REFERENCE;
    PREPARE_PATCHED_REFERENCE;
    L2_ALIGN;
    BUILD_NT_COUNTS;
    MERGE_SUMMARY_TABLE as MERGE_PREQC_SUMMARY;
    MERGE_SUMMARY_TABLE as MERGE_POSTQF_SUMMARY;
    MERGE_SUMMARY_TABLE as MERGE_POSTSUB_SUMMARY;
    MERGE_SUMMARY_TABLE as MERGE_POSTAT_SUMMARY;
    MERGE_SUMMARY_TABLE as MERGE_NONHUMAN_SUMMARY;
    MERGE_SUMMARY_TABLE as MERGE_BEST_REF_SUMMARY;
    MERGE_SUMMARY_TABLE as MERGE_Q40_SUMMARY;
    MERGE_SUMMARY_TABLE as MERGE_TARGET_SUMMARY;
    MERGE_SUMMARY_TABLE as MERGE_ASSEMBLY_READ_SUMMARY;
    MERGE_SUMMARY_TABLE as MERGE_CONTIG_COUNT_SUMMARY;
    MERGE_SUMMARY_TABLE as MERGE_L1_FLAG_SUMMARY;
    MERGE_SUMMARY_TABLE as MERGE_L2_FLAG_SUMMARY;
    MERGE_SUMMARY_TABLE as MERGE_L1_PRIMARY_CIGAR_SUMMARY;
    MERGE_SUMMARY_TABLE as MERGE_L2_PRIMARY_CIGAR_SUMMARY;
    MERGE_SUMMARY_TABLE as MERGE_L1_COVERAGE_SUMMARY;
    MERGE_SUMMARY_TABLE as MERGE_L2_COVERAGE_SUMMARY;
    MERGE_SUMMARY_TABLE as MERGE_L1_VARIANT_SUMMARY;
    MERGE_SUMMARY_TABLE as MERGE_L2_VARIANT_SUMMARY;
    MERGE_SUMMARY_TABLE as MERGE_L1_SAMTOOLSSEQ_SUMMARY;
    MERGE_SUMMARY_TABLE as MERGE_L1_BCFSEQ_SUMMARY;
    MERGE_SUMMARY_TABLE as MERGE_L2_SAMTOOLSSEQ_SUMMARY;
    MERGE_FASTA_COLLECTION as MERGE_L1_CALLEDBASES_FASTA;
    MERGE_FASTA_COLLECTION as MERGE_L1_BCF_FASTA;
    MERGE_FASTA_COLLECTION as MERGE_L2_CALLEDBASES_FASTA;
    MERGE_FASTA_COLLECTION as MERGE_L2_CONSENSUS_FASTA;
    MERGE_FASTA_COLLECTION as MERGE_L2_BCF_FASTA
} from '../modules/local/mlps_pe_processes'

workflow MLPS_PE_PIPELINE {
    take:
    input_samples
    primer_files

    main:
    project_name = new File(params.outdir.toString()).name

    // Step 1: normalize raw input names and expose canonical paired FASTQ filenames.
    // Output: `linked_fastqs.pair` -> tuple(meta, canonical_r1, canonical_r2)
    linked_fastqs = NORMALIZE_INPUT_FASTQ(input_samples)

    // Step 2: produce pre-QC read counts and visual QC reports from the canonical FASTQ files.
    // Output: summary TSV files plus per-sample PreQC TSV/PDF reports.
    preqc_counts = COUNT_PREQC_READS(
        linked_fastqs.pair.map { meta, read1, read2 -> tuple(meta, read1, read2, '00_preQC_read_counts') }
    )
    preqc_reports = GENERATE_PREQC(linked_fastqs.pair)

    // Step 3: quality filter, subsample, and adapter trim to prepare reads for mapping.
    // Output: `adapter_trimmed_fastqs.pair` -> cleaned per-sample paired FASTQ files.
    trimmed_fastqs = QUALITY_TRIM(linked_fastqs.pair)
    postqf_counts = COUNT_POSTQF_READS(
        trimmed_fastqs.pair.map { meta, read1, read2 -> tuple(meta, read1, read2, '02_postQF_read_counts') }
    )
    subsampled_fastqs = SUBSAMPLE_READS(trimmed_fastqs.pair)
    postsub_counts = COUNT_POSTSUB_READS(
        subsampled_fastqs.pair.map { meta, read1, read2 -> tuple(meta, read1, read2, '02_postSUB_read_counts') }
    )
    adapter_trimmed_fastqs = ADAPTER_TRIM(
        subsampled_fastqs.pair
            .combine(primer_files)
            .map { meta, read1, read2, fwdPrimer, revPrimer -> tuple(meta, read1, read2, fwdPrimer, revPrimer) }
    )
    postat_counts = COUNT_POSTAT_READS(
        adapter_trimmed_fastqs.pair.map { meta, read1, read2 -> tuple(meta, read1, read2, '03_postAT_read_counts') }
    )

    // Step 4: optionally remove human reads before target analysis.
    // Output: `nonhuman_fastqs` -> host-depleted paired FASTQ files when enabled.
    if( params.remove_human ) {
        human_reference_ready = CHECK_HUMAN_REFERENCE(Channel.value(params.human_db.toString()))
        nonhuman_fastqs = REMOVE_HUMAN_READS(
            adapter_trimmed_fastqs.pair
                .combine(human_reference_ready.check)
                .map { meta, read1, read2, ignored -> tuple(meta, read1, read2) }
        ).pair
    }
    else {
        nonhuman_fastqs = adapter_trimmed_fastqs.pair
    }

    // Step 5: capture post-filter read counts for downstream summaries.
    // Output: `nonhuman_counts.count_file` -> per-sample count TSV files.
    nonhuman_counts = COUNT_NONHUMAN_READS(
        nonhuman_fastqs.map { meta, read1, read2 -> tuple(meta, read1, read2, '05_nonhuman_counts') }
    )

    // Step 6: prepare formatted and doubled target references for linear/circular mapping.
    // Output: formatted and doubled FASTA assets.
    reference_assets = FORMAT_TARGET_REFERENCE(Channel.value(file(params.target_reference)))

    // Step 7: align reads to the mixed reference set, then estimate the dominant reference.
    // Output: mixed-reference BAMs, EM assignments, and per-sample selected FASTA references.
    mixed_bam_input = nonhuman_fastqs
        .combine(reference_assets.refs)
        .map { meta, read1, read2, refFormat, refDouble -> tuple(meta, read1, read2, refFormat, refDouble) }

    mixed_bams = MIXGENOME_ALIGN(mixed_bam_input)
    em_outputs = EM_REALLOCATE_REFERENCES(mixed_bams.bam)
    top_reference_input = em_outputs.summary
        .combine(reference_assets.formatted)
        .map { meta, summaryFile, formattedRef -> tuple(meta, summaryFile, formattedRef) }
    top_references = EXTRACT_TOP_REFERENCE(top_reference_input)
    sample_references = PREPARE_SAMPLE_REFERENCE(top_references.reference)
    sample_reference_by_id = sample_references.refs.map { meta, sampleRef, sampleDoubleRef -> tuple(meta.id, meta, sampleRef, sampleDoubleRef) }
    best_ref_rows = WRITE_BEST_REFERENCE(top_references.reference)

    // Step 8: extract target-mapping reads and build the L1 alignment set.
    // Output: target-only paired FASTQ files, L1 BAMs, and alignment flag summaries.
    target_read_pool = BUILD_TARGET_READ_POOL(mixed_bams.bam)
    q40_counts = COUNT_Q40_PRIMARY_READS(mixed_bams.bam)
    target_counts = COUNT_TARGET_READS(
        target_read_pool.pair.map { meta, read1, read2 -> tuple(meta, read1, read2, '06_potential_read_counts') }
    )
    target_reads_by_id = target_read_pool.pair.map { meta, read1, read2 -> tuple(meta.id, meta, read1, read2) }

    l1_align_input = target_reads_by_id
        .join(sample_reference_by_id)
        .map { sampleId, readMeta, read1, read2, refMeta, sampleRef, sampleDoubleRef -> tuple(readMeta, read1, read2, sampleRef, sampleDoubleRef) }
    l1_bams = L1_ALIGN(l1_align_input)
    l1_bams_by_id = l1_bams.bam.map { meta, bam -> tuple(meta.id, meta, bam) }
    l1_flag_counts = COUNT_L1_ALIGNMENT_FLAGS(l1_bams.bam)
    l1_primary_cigar = SUMMARISE_L1_PRIMARY_CIGAR(
        l1_bams.bam.map { meta, bam -> tuple(meta, bam, '07_L1_results', 'L1') }
    )
    l1_depth_input = l1_bams_by_id
        .join(sample_reference_by_id)
        .map { sampleId, bamMeta, bam, refMeta, sampleRef, sampleDoubleRef -> tuple(bamMeta, bam, sampleRef, '07_L1_results') }
    l1_depth_mask = BUILD_L1_DEPTH_MASK(l1_depth_input)
    l1_coverage_rows = SUMMARISE_L1_COVERAGE(l1_depth_mask.depth_mask.map { meta, depthFile, maskFile -> tuple(meta, depthFile, '07_L1_results') })
    BUILD_L1_COVERAGE_PLOT(
        l1_depth_mask.depth_mask.map { meta, depthFile, maskFile -> depthFile }.collect(),
        "${project_name}.L1dedup.coverage.html"
    )

    // Step 9: convert L1 alignments into pileup-driven consensus outputs.
    // Output: nt-frequency CSVs plus L1 consensus/reference FASTA files.
    l1_pileup_input = l1_bams_by_id
        .join(sample_reference_by_id)
        .map { sampleId, bamMeta, bam, refMeta, sampleRef, sampleDoubleRef -> tuple(bamMeta, bam, sampleRef, sampleDoubleRef) }
    l1_pileups = BUILD_L1_PILEUP(l1_pileup_input)
    l1_consensus = FORMAT_L1_PILEUP(
        l1_pileups.pileup.map { meta, pileup, original_ref -> tuple(meta, pileup, original_ref, '07_L1_results') }
    )
    l1_variant_input = l1_bams_by_id
        .join(sample_reference_by_id)
        .join(l1_depth_mask.depth_mask.map { meta, depthFile, maskFile -> tuple(meta.id, meta, depthFile, maskFile) })
        .map { sampleId, bamMeta, bam, refMeta, sampleRef, sampleDoubleRef, maskMeta, depthFile, maskFile -> tuple(bamMeta, bam, sampleRef, sampleDoubleRef, maskFile, '07_L1_results') }
    l1_variants = BUILD_L1_VARIANTS(l1_variant_input)

    // Step 10: assemble contigs from normalized target reads and project them back to the reference.
    // Output: assembly contigs, assembly BAMs, and assembly consensus FASTA files.
    normalized_for_assembly = NORMALISE_FOR_ASSEMBLY(target_read_pool.pair)
    assembly_read_counts = COUNT_ASSEMBLY_READS(
        normalized_for_assembly.pair.map { meta, read1, read2 -> tuple(meta, read1, read2, '08_assembly_read_counts') }
    )
    spades_contigs = RUN_SPADES(normalized_for_assembly.pair)
    contig_counts = COUNT_CONTIG_READS(
        spades_contigs.contigs.map { meta, contigs -> tuple(meta, contigs, '08_contig_counts') }
    )
    spades_contigs_by_id = spades_contigs.contigs.map { meta, contigs -> tuple(meta.id, meta, contigs) }
    assembly_align_input = spades_contigs_by_id
        .join(sample_reference_by_id)
        .map { sampleId, contigMeta, contigs, refMeta, sampleRef, sampleDoubleRef -> tuple(contigMeta, contigs, sampleRef, sampleDoubleRef) }
    assembly_bams = ALIGN_ASSEMBLY_TO_REFERENCE(assembly_align_input)
    assembly_bams_by_id = assembly_bams.bam.map { meta, bam -> tuple(meta.id, meta, bam) }
    assembly_pileup_input = assembly_bams_by_id
        .join(sample_reference_by_id)
        .map { sampleId, bamMeta, bam, refMeta, sampleRef, sampleDoubleRef -> tuple(bamMeta, bam, sampleRef, sampleDoubleRef) }
    assembly_pileups = BUILD_ASSEMBLY_PILEUP(assembly_pileup_input)
    assembly_consensus = FORMAT_ASSEMBLY_PILEUP(
        assembly_pileups.pileup.map { meta, pileup, original_ref -> tuple(meta, pileup, original_ref, '08_assembly') }
    )

    // Step 11: patch the sample reference with assembly evidence and realign for L2 consensus.
    // Output: patched reference FASTA files plus final L2 BAMs and consensus FASTA files.
    patched_references = PREPARE_PATCHED_REFERENCE(assembly_consensus.patched_reference)
    patched_reference_by_id = patched_references.refs.map { meta, originalRef, mappingRef -> tuple(meta.id, meta, originalRef, mappingRef) }
    l2_align_input = target_reads_by_id
        .join(patched_reference_by_id)
        .map { sampleId, readMeta, read1, read2, refMeta, originalRef, mappingRef -> tuple(readMeta, read1, read2, originalRef, mappingRef) }
    l2_bams = L2_ALIGN(l2_align_input)
    l2_bams_by_id = l2_bams.bam.map { meta, bam -> tuple(meta.id, meta, bam) }
    l2_flag_counts = COUNT_L2_ALIGNMENT_FLAGS(l2_bams.bam)
    l2_primary_cigar = SUMMARISE_L2_PRIMARY_CIGAR(
        l2_bams.bam.map { meta, bam -> tuple(meta, bam, '09_L2_results', 'L2') }
    )
    l2_depth_input = l2_bams_by_id
        .join(patched_reference_by_id)
        .map { sampleId, bamMeta, bam, refMeta, originalRef, mappingRef -> tuple(bamMeta, bam, originalRef, '09_L2_results') }
    l2_depth_mask = BUILD_L2_DEPTH_MASK(l2_depth_input)
    l2_coverage_rows = SUMMARISE_L2_COVERAGE(l2_depth_mask.depth_mask.map { meta, depthFile, maskFile -> tuple(meta, depthFile, '09_L2_results') })
    BUILD_L2_COVERAGE_PLOT(
        l2_depth_mask.depth_mask.map { meta, depthFile, maskFile -> depthFile }.collect(),
        "${project_name}.L2dedup.coverage.html"
    )
    l2_pileup_input = l2_bams_by_id
        .join(patched_reference_by_id)
        .map { sampleId, bamMeta, bam, refMeta, originalRef, mappingRef -> tuple(bamMeta, bam, originalRef, mappingRef) }
    l2_pileups = BUILD_L2_PILEUP(l2_pileup_input)
    l2_consensus = FORMAT_L2_PILEUP(
        l2_pileups.pileup.map { meta, pileup, original_ref -> tuple(meta, pileup, original_ref, '09_L2_results') }
    )
    if( params.target_genome_type == 'linear' ) {
        l2_variant_input = l2_bams_by_id
            .join(patched_reference_by_id)
            .join(l2_depth_mask.depth_mask.map { meta, depthFile, maskFile -> tuple(meta.id, meta, depthFile, maskFile) })
            .map { sampleId, bamMeta, bam, refMeta, originalRef, mappingRef, maskMeta, depthFile, maskFile -> tuple(bamMeta, bam, originalRef, mappingRef, maskFile, '09_L2_results') }
        l2_variant_rows = BUILD_L2_VARIANTS(l2_variant_input)
    }
    else {
        l2_variant_rows = WRITE_VARIANT_STATUS_PLACEHOLDER(
            l2_bams.bam.map { meta, bam -> tuple(meta, '09_L2_results') }
        )
    }

    // Step 12: expand non-aligned nucleotide summaries and merge the final report collections.
    // Output: per-position nt counts plus merged summary/FASTA deliverables.
    nt_counts = BUILD_NT_COUNTS(l2_consensus.ntfreq)

    MERGE_PREQC_SUMMARY(
        preqc_counts.count_file.map { meta, stage, countFile -> countFile }.collect(),
        '00_preQC_read_counts.tsv',
        'job_id\tPreQC_reads\tPreQC_bases'
    )
    MERGE_POSTQF_SUMMARY(
        postqf_counts.count_file.map { meta, stage, countFile -> countFile }.collect(),
        '02_postQF_read_counts.tsv',
        'job_id\tPostQF_reads\tPostQF_bases'
    )
    MERGE_POSTSUB_SUMMARY(
        postsub_counts.count_file.map { meta, stage, countFile -> countFile }.collect(),
        '02_postSUB_read_counts.tsv',
        'job_id\tPostSUB_reads\tPostSUB_bases'
    )
    MERGE_POSTAT_SUMMARY(
        postat_counts.count_file.map { meta, stage, countFile -> countFile }.collect(),
        '03_postAT_read_counts.tsv',
        'job_id\tPostAT_reads\tPostAT_bases'
    )
    MERGE_NONHUMAN_SUMMARY(
        nonhuman_counts.count_file.map { meta, stage, countFile -> countFile }.collect(),
        '05_nonhuman_counts.tsv',
        'job_id\tnonhuman_reads\tnonhuman_bases'
    )
    MERGE_BEST_REF_SUMMARY(
        best_ref_rows.summary_file.map { meta, summaryFile -> summaryFile }.collect(),
        '06_best_ref.tsv',
        'job_id\tBest_ref'
    )
    MERGE_Q40_SUMMARY(
        q40_counts.count_file.map { meta, countFile -> countFile }.collect(),
        '06_multi_q40_read_counts.tsv',
        'job_id\tMulti_q40_reads'
    )
    MERGE_TARGET_SUMMARY(
        target_counts.count_file.map { meta, stage, countFile -> countFile }.collect(),
        '06_potential_read_counts.tsv',
        'job_id\tPotential_reads\tPotential_bases'
    )
    MERGE_ASSEMBLY_READ_SUMMARY(
        assembly_read_counts.count_file.map { meta, stage, countFile -> countFile }.collect(),
        '08_assembly_read_counts.tsv',
        'job_id\tassembly_reads\tassembly_bases'
    )
    MERGE_CONTIG_COUNT_SUMMARY(
        contig_counts.count_file.map { meta, stage, countFile -> countFile }.collect(),
        '08_contig_counts.tsv',
        'job_id\tcontig_reads\tcontig_bases'
    )
    MERGE_L1_FLAG_SUMMARY(
        l1_flag_counts.flag_file.map { meta, flagFile -> flagFile }.collect(),
        '07_dedup_counts.tsv',
        'job_id\tL1_fwd\tL1_rev\tL1_fwd_dup\tL1_rev_dup'
    )
    MERGE_L2_FLAG_SUMMARY(
        l2_flag_counts.flag_file.map { meta, flagFile -> flagFile }.collect(),
        '09_dedup_counts.tsv',
        'job_id\tL2_fwd\tL2_rev\tL2_fwd_dup\tL2_rev_dup'
    )
    MERGE_L1_PRIMARY_CIGAR_SUMMARY(
        l1_primary_cigar.cigar_row.map { meta, rowFile -> rowFile }.collect(),
        '07_primarycigar.tsv',
        'job_id\tL1_hardclip_bases\tL1_softclip_bases\tL1_match_bases\tL1_insertion_bases\tL1_deletion_bases'
    )
    MERGE_L2_PRIMARY_CIGAR_SUMMARY(
        l2_primary_cigar.cigar_row.map { meta, rowFile -> rowFile }.collect(),
        '09_primarycigar.tsv',
        'job_id\tL2_hardclip_bases\tL2_softclip_bases\tL2_match_bases\tL2_insertion_bases\tL2_deletion_bases'
    )
    MERGE_L1_COVERAGE_SUMMARY(
        l1_coverage_rows.coverage_row.map { meta, rowFile -> rowFile }.collect(),
        '07_L1_coverage.tsv',
        'job_id\tave_covs\tmedian_covs\thighest_ave_cov\thighest_median_cov'
    )
    MERGE_L2_COVERAGE_SUMMARY(
        l2_coverage_rows.coverage_row.map { meta, rowFile -> rowFile }.collect(),
        '09_L2_coverage.tsv',
        'job_id\tave_covs\tmedian_covs\thighest_ave_cov\thighest_median_cov'
    )
    MERGE_L1_VARIANT_SUMMARY(
        l1_variants.status_file.map { meta, rowFile -> rowFile }.collect(),
        '07_variant_status.tsv',
        'job_id\tL1_variant_sites\tL1_high_qual_variants'
    )
    MERGE_L2_VARIANT_SUMMARY(
        l2_variant_rows.status_file.map { meta, rowFile -> rowFile }.collect(),
        '09_variant_status.tsv',
        'job_id\tL2_variant_sites\tL2_high_qual_variants'
    )
    MERGE_L1_SAMTOOLSSEQ_SUMMARY(
        l1_consensus.fasta_summary_row.map { meta, rowFile -> rowFile }.collect(),
        '07_samtoolsseq_summmary.tsv',
        'job_id\tL1_high_cof_pos\tL1_low_cof_pos\tL1_ambiguous_pos\tL1_gap_pos'
    )
    MERGE_L1_BCFSEQ_SUMMARY(
        l1_variants.fasta_summary_row.map { meta, rowFile -> rowFile }.collect(),
        '07_bcfseq_summmary.tsv',
        'job_id\tL1_high_cof_pos\tL1_low_cof_pos\tL1_ambiguous_pos\tL1_gap_pos'
    )
    MERGE_L2_SAMTOOLSSEQ_SUMMARY(
        l2_consensus.fasta_summary_row.map { meta, rowFile -> rowFile }.collect(),
        '09_samtoolsseq_summmary.tsv',
        'job_id\tL2_high_cof_pos\tL2_low_cof_pos\tL2_ambiguous_pos\tL2_gap_pos'
    )

    MERGE_L1_CALLEDBASES_FASTA(
        l1_consensus.bases_fasta.map { meta, fastaFile -> fastaFile }.collect(),
        "${project_name}.L1.calledbases.fasta"
    )
    MERGE_L1_BCF_FASTA(
        l1_variants.consensus_fasta.map { meta, fastaFile -> fastaFile }.collect(),
        "${project_name}.L1.bcfconsensus.fasta"
    )
    MERGE_L2_CALLEDBASES_FASTA(
        l2_consensus.bases_fasta.map { meta, fastaFile -> fastaFile }.collect(),
        "${project_name}.L2.calledbases.fasta"
    )
    MERGE_L2_CONSENSUS_FASTA(
        l2_consensus.consensus_fasta.map { meta, fastaFile -> fastaFile }.collect(),
        "${project_name}.L2.consensus.fasta"
    )
    if( params.target_genome_type == 'linear' ) {
        MERGE_L2_BCF_FASTA(
            l2_variant_rows.consensus_fasta.map { meta, fastaFile -> fastaFile }.collect(),
            "${project_name}.L2.bcfconsensus.fasta"
        )
    }

    emit:
    // Key emitted channels for debugging and reuse in future subworkflows.
    linked = linked_fastqs.pair
    nonhuman = nonhuman_fastqs
    l2 = l2_consensus.bases_fasta
    nt = nt_counts.nt_csv
}
