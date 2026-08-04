# MLPS_PE

`MLPS_PE` is a Nextflow pipeline for running the paired-end MLPS analysis workflow in a reproducible, production-ready form. It translates the existing shell-based process into a structured workflow that is easier to run, validate, and maintain across datasets and execution environments.

The pipeline is designed for short-read paired-end sequencing data and preserves the established MLPS analysis flow, from input preparation and quality filtering through host removal, target alignment, assembly, and result summarization. Its goal is to provide a consistent interface for routine analysis while keeping the expected result structure familiar to existing users of the workflow.

Project author: Dr. Haiting Chai

## Prerequisite

Install Nextflow before running the pipeline:

```bash
curl -s https://get.nextflow.io | bash
```

## MLPS_PE Command-Line Toolkit

The repository includes a terminal dispatcher for the Python toolkit. It uses `python/toolkit/toolkit.list.csv` as the module registry and provides one command for running the registered utilities.

Show the toolkit overview and all available modules:

```bash
./mlps_pe --help
```

Show detailed help for one module. The module options are forwarded to the underlying Python script:

```bash
./mlps_pe format_fasta --help
```

Run a module:

```bash
./mlps_pe format_fasta \
  --input input.fa \
  --output output.fa
```

The dispatcher and toolkit registry are located in `python/toolkit/`:

- `python/toolkit/mlps_pe.py`: module dispatcher
- `python/toolkit/toolkit.list.csv`: registered modules and parameter summaries
- `python/toolkit/requirements-toolkit.txt`: Python package requirements

### Building the Toolkit Runtime

On macOS or Linux, run the build script from the project root:

```bash
./make
```

This creates `.mlps_pe_runtime/` and installs the toolkit’s Python package dependencies. The existing `./mlps_pe` launcher automatically uses that runtime.

On Windows PowerShell, run:

```powershell
.\make.ps1
```

The Windows build creates `.mlps_pe_runtime\` and installs the same package requirements. The dispatcher can then be run with:

```powershell
.\.mlps_pe_runtime\Scripts\python.exe .\python\toolkit\mlps_pe.py --help
```

The `.mlps_pe_runtime/` directory is a local generated runtime and is excluded from version control.

## Workflow Scope

This first-pass production pipeline follows the process order documented in `Pilot/workflow.md` and preserves the recognizable numbered result layout from the shell workflow:

- `00_fastq_links`
- `01_fastq_reports`
- `02_quality_trimming`
- `03_adapter_trimming`
- `04_reference`
- `05_nonhuman`
- `06_mixgenome`
- `07_L1_results`
- `08_assembly`
- `09_L2_results`
- `10_nonaligned_nt_counts`
- `summary`

The pipeline reuses the repository’s helper scripts in `python/` and the conda YAMLs in `environment/`.

## Supported Inputs

- paired-end short reads only
- gzipped FASTQ input paths supplied through `sample_list.csv`
- short reads only, currently expected to be `<500 bp`

The current production pipeline is designed for short-read paired-end data and is not validated for reads `>=500 bp`.

The input CSV must contain these header columns:

- `job_id`
- `sample_name`
- `r1_path`
- `r2_path`

Relative `r1_path` and `r2_path` values are resolved relative to the CSV file location.

Canonical internal naming now follows the paired-end workflow design:

- R1 is linked as `${job_id}.R1.fastq.gz`
- R2 is linked as `${job_id}.R2.fastq.gz`

## Required Parameters

The pipeline validates parameters and key repository assets before expensive work begins.

- `--inputcsv`: input sample sheet CSV
- `--outdir`: output directory
- `--human_db`: Bowtie2 human index prefix
- `--target_reference`: target reference FASTA
- `--fwd_primer_seq`: forward primer sequence to trim in R1 file
- `--rev_primer_seq`: reverse primer sequence to trim in R2 file
- `--dedup`: `yes` or `no`
- `--target_genome_type`: `linear` or `circular`
- `--bq_thr`: base-quality threshold for bcftools mpileup, default `13`
- `--depth_thr`: depth threshold used for masking and filtering, default `5`
- `--qual_thr`: QUAL threshold for filtered variants, default `30`
- `--umi_location`:location of the unique molecular identifier, `read1`, `read2` or `per_read`, default `none`
- `--umi_length`: length of the UMI, default `8`
- `--min_insert`: minimum post-QC length of the insert, default `30`
- `--cleantmpdata`: `yes` or `no`, default `yes`
- `--cleantmpdatasize`: size threshold for truncating oversized regular files in traced `work/` task directories during shutdown cleanup, default `1MB`
- `--threads`: positive integer

The human-index example below is local-only documentation and is not committed as a default:

```bash
/Users/macbook/Documents/HumanGenome/T2T/t2t_index
```

## Example Usage

Example `sample_list.csv`:

```csv
job_id,sample_name,r1_path,r2_path
test001,HBV1,./test/linear/HBV1.R1.fastq.gz,./test/linear/HBV1.R2.fastq.gz
test002,HBV2,./test/linear/HBV2.R1.fastq.gz,./test/linear/HBV2.R2.fastq.gz
```

Linear example:

```bash
nextflow run . \
  --inputcsv ./sample_list.csv \
  --outdir ./results/hbv_linear_run \
  --human_db <human_db_path>/<db_prefix> \
  --target_reference ./reference/ref.HBV.04092023.fasta \
  --target_genome_type linear \
  --dedup yes \
  --bq_thr 13 \
  --depth_thr 5 \
  --qual_thr 30 \
  --threads 4
```

Circular example:

```bash
nextflow run . \
  --inputcsv ./sample_list.csv \
  --outdir ./results/hbv_run \
  --human_db <human_db_path>/<db_prefix> \
  --target_reference ./reference/ref.HBV.04092023.fasta \
  --target_genome_type circular \
  --dedup yes \
  --bq_thr 13 \
  --depth_thr 5 \
  --qual_thr 30 \
  --threads 4
```

## Profiles

- `local`: workstation execution
- `test`: lightweight linear smoke defaults using bundled HBV test data
- `ci`: CI-friendly settings that disable conda solving during stub runs

## Bundled Test Data

- `test/linear/HBV1.R1.fastq.gz`
- `test/linear/HBV1.R2.fastq.gz`
- `test/linear/HBV2.R1.fastq.gz`
- `test/linear/HBV2.R2.fastq.gz`
- `test/circular/HBV1.R1.fastq.gz`
- `test/circular/HBV1.R2.fastq.gz`
- `sample_list.csv`

Linear smoke test:

```bash
nextflow run . -profile test,ci -stub-run -ansi-log false
```

Circular smoke test:

```bash
nextflow run . \
  -profile test,ci \
  -stub-run \
  -ansi-log false \
  --inputcsv ./sample_list.csv \
  --outdir ./results/test_circular \
  --target_reference ./reference/ref.HBV.04092023.fasta \
  --target_genome_type circular
```

First-time setup recommendation:

- Run the bundled test samples first with the small bundled human database so Nextflow can create and cache the required conda environments before you process real data.
- After those environments have been created successfully on a machine with internet access, later runs can usually reuse the cached environments in offline mode as long as the same environment definitions and cache locations are available.
- `-stub-run` validates workflow structure only. To actually populate the conda environment cache, run the test configuration without `-stub-run`.

Example first real test run for environment caching:

```bash
nextflow run . -c conf/test.config --remove_human yes
```

## Key Published Outputs

Important documented summary outputs include:

- `01_fastq_reports/PreQC_PDF/${sample_name}.preQC.pdf`
- `06_mixgenome/${sample_name}.em_assignments_ref_summary.tsv`
- `07_L1_results/${sample_name}.primarycigar.txt`
- `07_L1_results/${sample_name}.norm.vcf.tsv`
- `09_L2_results/${sample_name}.primarycigar.txt`
- `09_L2_results/${sample_name}.norm.vcf.tsv` for linear runs
- `10_nonaligned_nt_counts/${sample_name}.L2nt.csv`
- `summary/${project_name}.results.tsv`, generated at workflow shutdown from all available summary TSVs so partial runs still emit a merged table
- `${project_name}.L1.calledbases.fasta`
- `${project_name}.L1.bcfconsensus.fasta`
- `${project_name}.L2.calledbases.fasta`
- `${project_name}.L2.consensus.fasta`
- `${project_name}.L2.bcfconsensus.fasta` for linear runs
- `pipeline_info/${project_name}.command.md`, generated at workflow shutdown from available Nextflow `.command.sh` task scripts for command tracing

When `--cleantmpdata yes` is active, the shutdown finalizer also reads `pipeline_info/trace.txt` and truncates oversized regular files inside the traced task directories under `work/`. It skips symlinks and does not touch regular files outside the Nextflow `work/` tree.

## CI

GitHub Actions under `.github/workflows/ci.yml` currently performs:

- repository checkout
- Python setup
- Java setup
- Python syntax checking, including `python/toolkit/mlps_pe.py`
- MLPS toolkit dispatcher help and registry smoke check via `./mlps_pe --help`
- Nextflow installation
- config validation
- paired-end linear stub smoke execution
- paired-end circular stub smoke execution with a temporary circular sample sheet

## Current Validation Status

This repository currently validates well in stub mode and preserves the documented production layout. Real non-stub biological execution still depends on solving and running the conda environments with the required external tools available.
