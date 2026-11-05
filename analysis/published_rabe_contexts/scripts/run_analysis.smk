import csv
import re

# Run from the repository root. Biological metadata belongs in the CSV inputs.
ANALYSIS = "analysis/published_rabe_contexts"
SCRIPTS = f"{ANALYSIS}/scripts"
ANNOTATIONS = f"{ANALYSIS}/annotations"
DATA = f"{ANALYSIS}/data"
TABLES = f"{ANALYSIS}/tables"
LOGS = f"{ANALYSIS}/logs"
BENCHMARKS = f"{ANALYSIS}/benchmarks"
SAMPLES_FILE = f"{ANNOTATIONS}/sample_annotations.csv"
RUNS_FILE = f"{ANNOTATIONS}/run_annotations.csv"
REFERENCE_FILE = f"{ANNOTATIONS}/reference_annotations.csv"
CONTAINERS_FILE = f"{ANNOTATIONS}/containers.csv"
READ_MODE = "read1"  # Read 1 only, of the paired-end libraries.
RNA_STRANDNESS = "F"
BASE_CHANGE = "A,G"
REPEAT_LIMIT = 1000
MIN_COVERAGE_EXCLUSIVE = 10
CHUNK_SIZE = 100000
SORT_MEMORY = "1G"
SORT_SCRATCH = f"{DATA}/sort_scratch"
ALIGN_THREADS = 12
SORT_EXTRA_THREADS = 3
TABLE_THREADS = 10
BUILD_THREADS = 16
BUILD_MEMORY = 300000
BUILD_RUNTIME = 1440
FASTA = f"{DATA}/reference/genome.fa"
GTF = f"{DATA}/reference/annotations.gtf"
INDEX = f"{DATA}/reference/index/genome"
INDEX_FILES = [f"{INDEX}.3n.{alphabet}{repeat}.{part}.ht2"
               for alphabet in ("CT", "GA")
               for repeat, parts in (("", range(1, 9)), (".rep", (1, 2, 3, 4, 7, 8)))
               for part in parts]
FIGURE_NAME = "rabe_rna_sequence_context"
FIGURE_PNG = f"{ANALYSIS}/figures/{FIGURE_NAME}.png"
# The manuscript panel this figure is, exported beside every other panel of the paper.
PANEL_CSV = "source_data/figure_s3b.csv"
FIGURE_PDF = f"{ANALYSIS}/figures/{FIGURE_NAME}.pdf"

if config:
    raise ValueError("Use the full workflow without configuration overrides; choices are declared above")

with open(SAMPLES_FILE, newline="") as handle:
    sample_rows = list(csv.DictReader(handle))
with open(RUNS_FILE, newline="") as handle:
    run_rows = list(csv.DictReader(handle))
with open(REFERENCE_FILE, newline="") as handle:
    reference_rows = list(csv.DictReader(handle))
with open(CONTAINERS_FILE, newline="") as handle:
    containers = {row["tool"]: row["uri"] for row in csv.DictReader(handle)}

SAMPLES = [row["sample_name"] for row in sample_rows]
RUNS = {row["run_accession"]: row for row in run_rows}
SAMPLE_RUNS = {sample: [row["run_accession"] for row in run_rows
                       if row["sample_name"] == sample] for sample in SAMPLES}
if not SAMPLES or len(set(SAMPLES)) != len(SAMPLES) or len(RUNS) != len(run_rows):
    raise ValueError("Sample names and run accessions must be nonempty and unique")
if len(reference_rows) != 1:
    raise ValueError("Exactly one reference annotation row is required")
if READ_MODE != "read1":
    raise ValueError("Only the approved read-1-only analysis is implemented")
if set(SAMPLES) != {row["sample_name"] for row in run_rows}:
    raise ValueError("Every annotated sample must have runs; unexpected samples are forbidden")
if any(not re.fullmatch(r"[A-Za-z0-9_]+", sample) for sample in SAMPLES):
    raise ValueError("Sample names must be safe filename components")
for sample in sample_rows:
    if not sample["condition"] or not sample["replicate"].isdigit():
        raise ValueError("Every sample needs an explicit condition and replicate")
    if not sample["condition_order"].isdigit():
        raise ValueError("Every condition needs a numeric display order")
    if any(RUNS[run]["experiment_accession"] != sample["experiment_accession"]
           for run in SAMPLE_RUNS[sample["sample_name"]]):
        raise ValueError("Run and sample experiment accessions disagree")
if len({(row["condition"], row["replicate"]) for row in sample_rows}) != len(SAMPLES):
    raise ValueError("Duplicate condition/replicate annotation")
condition_orders = {(row["condition"], row["condition_order"]) for row in sample_rows}
if len(condition_orders) != len({row["condition"] for row in sample_rows}):
    raise ValueError("Inconsistent display order within a condition")
if len(condition_orders) != len({row["condition_order"] for row in sample_rows}):
    raise ValueError("Different conditions share a display order")
for row in run_rows:
    if row["library_layout"] != "PAIRED":
        raise ValueError("Unexpected source layout: review read selection before running")
    if not re.fullmatch(r"[SED]RR[0-9]+", row["run_accession"]):
        raise ValueError("Invalid run accession")
    for mate in (1, 2):
        if not row[f"read{mate}_url"].startswith("https://ftp.sra.ebi.ac.uk/"):
            raise ValueError("Unexpected FASTQ download host")
        if not re.fullmatch(r"[0-9a-f]{32}", row[f"read{mate}_md5"]):
            raise ValueError("Missing FASTQ checksum")
REFERENCE = reference_rows[0]

wildcard_constraints:
    sample="|".join(re.escape(sample) for sample in SAMPLES),
    run="|".join(re.escape(run) for run in RUNS),
    mate="[12]"


rule all:
    input:
        FIGURE_PNG,
        FIGURE_PDF,
        f"{TABLES}/{FIGURE_NAME}.csv",
        PANEL_CSV,
        expand(f"{TABLES}/qc/{{sample}}.csv", sample=SAMPLES),
        expand(f"{DATA}/sites/{{sample}}.csv.gz", sample=SAMPLES)


rule download_fastq:
    input:
        RUNS_FILE
    output:
        fastq=f"{DATA}/fastq/{{run}}_{{mate}}.fastq.gz"
    params:
        url=lambda w: RUNS[w.run][f"read{w.mate}_url"],
        checksum=lambda w: RUNS[w.run][f"read{w.mate}_md5"]
    threads: 1
    resources:
        mem_mb=1000, runtime=240, download_slots=1
    container: containers["analysis"]
    log: f"{LOGS}/download/{{run}}_{{mate}}.log"
    shell:
        r"""
        curl --fail --location --retry 3 --retry-delay 10 --output {output.fastq:q}.partial \
            {params.url:q} > {log:q} 2>&1
        printf '%s  %s\n' {params.checksum:q} {output.fastq:q}.partial | md5sum --check >> {log:q} 2>&1
        gzip -t {output.fastq:q}.partial 2>> {log:q}
        mv {output.fastq:q}.partial {output.fastq:q}
        """


rule download_reference:
    input:
        REFERENCE_FILE
    output:
        fasta=FASTA, gtf=GTF
    params:
        fasta_url=REFERENCE["fasta_url"], gtf_url=REFERENCE["gtf_url"]
    threads: 1
    resources:
        mem_mb=1000, runtime=240, download_slots=1
    container: containers["analysis"]
    log: f"{LOGS}/reference/download.log"
    shell:
        r"""
        curl --fail --location --retry 3 {params.fasta_url:q} 2> {log:q} | gzip -dc > {output.fasta:q}
        curl --fail --location --retry 3 {params.gtf_url:q} 2>> {log:q} | gzip -dc > {output.gtf:q}
        test -s {output.fasta:q}
        test -s {output.gtf:q}
        """


rule index_fasta:
    input:
        fasta=FASTA
    output:
        fai=f"{FASTA}.fai"
    threads: 1
    resources:
        mem_mb=2000, runtime=30
    container: containers["alignment"]
    log: f"{LOGS}/reference/faidx.log"
    shell:
        "samtools faidx {input.fasta:q} > {log:q} 2>&1"


rule splice_annotations:
    input:
        gtf=GTF
    output:
        exons=f"{DATA}/reference/exons.tsv",
        splice_sites=f"{DATA}/reference/splice_sites.tsv"
    threads: 1
    resources:
        mem_mb=8000, runtime=60
    container: containers["alignment"]
    log: f"{LOGS}/reference/splicing.log"
    shell:
        r"""
        extract_exons.py {input.gtf:q} > {output.exons:q} 2> {log:q}
        extract_splice_sites.py {input.gtf:q} > {output.splice_sites:q} 2>> {log:q}
        """


rule build_index:
    input:
        fasta=FASTA,
        exons=rules.splice_annotations.output.exons,
        splice_sites=rules.splice_annotations.output.splice_sites
    output:
        INDEX_FILES
    params:
        index=lambda w, output: output[0].split(".3n.")[0], base_change=BASE_CHANGE
    threads: BUILD_THREADS
    resources:
        mem_mb=BUILD_MEMORY, runtime=BUILD_RUNTIME
    container: containers["alignment"]
    log: f"{LOGS}/reference/build_index.log"
    benchmark: f"{BENCHMARKS}/build_index.tsv"
    shell:
        r"""
        hisat-3n-build -p {threads} --base-change {params.base_change:q} --repeat-index \
            --exon {input.exons:q} --ss {input.splice_sites:q} \
            {input.fasta:q} {params.index:q} > {log:q} 2>&1
        for file in {output:q}; do
            test -s "$file" || {{ echo "Missing or empty index component: $file" >> {log:q}; exit 1; }}
        done
        """


rule align_genome:
    input:
        fastqs=lambda w: [f"{DATA}/fastq/{run}_1.fastq.gz" for run in SAMPLE_RUNS[w.sample]],
        indices=INDEX_FILES
    output:
        bam=f"{DATA}/alignments/{{sample}}.bam",
        bai=f"{DATA}/alignments/{{sample}}.bam.bai"
    params:
        index=lambda w, input: input.indices[0].split(".3n.")[0],
        base_change=BASE_CHANGE, repeat_limit=REPEAT_LIMIT,
        strandness=RNA_STRANDNESS,
        reads=lambda w, input: ",".join(input.fastqs),
        align_threads=lambda w, threads: max(1, threads - SORT_EXTRA_THREADS - 1),
        sort_threads=lambda w, threads: min(SORT_EXTRA_THREADS, max(0, threads - 2)),
        sort_memory=SORT_MEMORY, scratch_root=SORT_SCRATCH
    threads: ALIGN_THREADS
    resources:
        mem_mb=32000, runtime=720
    container: containers["alignment"]
    log: f"{LOGS}/alignments/{{sample}}.log"
    benchmark: f"{BENCHMARKS}/alignments/{{sample}}.tsv"
    shell:
        r"""
        mkdir -p {params.scratch_root:q}
        scratch=$(mktemp -d {params.scratch_root:q}/{wildcards.sample:q}.XXXXXX)
        trap 'rm -rf "$scratch"' EXIT
        hisat-3n -x {params.index:q} -p {params.align_threads} -q -U {params.reads:q} \
            --repeat --repeat-limit {params.repeat_limit} --base-change {params.base_change:q} \
            --no-unal --rna-strandness {params.strandness:q} 2> {log:q} | \
            samtools sort -@ {params.sort_threads} -m {params.sort_memory:q} \
                -T "$scratch/sort" -o {output.bam:q} - 2>> {log:q}
        samtools quickcheck -v {output.bam:q} >> {log:q} 2>&1
        samtools index {output.bam:q} {output.bai:q} >> {log:q} 2>&1
        """


rule count_conversions:
    input:
        bam=rules.align_genome.output.bam,
        fasta=FASTA,
        script=f"{SCRIPTS}/compact_counts.py"
    output:
        counts=f"{DATA}/counts/{{sample}}.csv.gz"
    params:
        base_change=BASE_CHANGE,
        table_threads=lambda w, threads: max(1, threads - 2)
    threads: TABLE_THREADS
    resources:
        mem_mb=32000, runtime=720
    container: containers["alignment"]
    log: f"{LOGS}/counts/{{sample}}.log"
    benchmark: f"{BENCHMARKS}/counts/{{sample}}.tsv"
    shell:
        r"""
        {{ samtools view -h {input.bam:q} | \
            hisat-3n-table --threads {params.table_threads} --alignments - \
                --ref {input.fasta:q} --base-change {params.base_change:q} | \
            python {input.script:q} {output.counts:q}; }} 2> {log:q}
        """


rule extract_context:
    input:
        counts=rules.count_conversions.output.counts,
        fasta=FASTA,
        fai=rules.index_fasta.output.fai,
        script=f"{SCRIPTS}/extract_sequence_context.r"
    output:
        sites=f"{DATA}/sites/{{sample}}.csv.gz",
        summary=f"{TABLES}/samples/{{sample}}.csv",
        qc=f"{TABLES}/qc/{{sample}}.csv"
    params:
        cutoff=MIN_COVERAGE_EXCLUSIVE, chunk_size=CHUNK_SIZE,
        converted_base=BASE_CHANGE.split(",")[0]
    threads: 1
    resources:
        mem_mb=8000, runtime=240
    container: containers["analysis"]
    log: f"{LOGS}/context/{{sample}}.log"
    benchmark: f"{BENCHMARKS}/context/{{sample}}.tsv"
    shell:
        r"""
        Rscript {input.script:q} {input.fasta:q} {input.counts:q} {output.sites:q} \
            {output.summary:q} {output.qc:q} {wildcards.sample:q} \
            {params.cutoff} {params.chunk_size} {params.converted_base:q} > {log:q} 2>&1
        """


rule plot_context:
    input:
        annotations=SAMPLES_FILE,
        summaries=expand(f"{TABLES}/samples/{{sample}}.csv", sample=SAMPLES),
        script=f"{SCRIPTS}/plot_context.r"
    output:
        source=f"{TABLES}/{FIGURE_NAME}.csv",
        panel=PANEL_CSV,
        png=FIGURE_PNG, pdf=FIGURE_PDF
    threads: 1
    resources:
        mem_mb=2000, runtime=30
    container: containers["analysis"]
    log: f"{LOGS}/plot_context.log"
    shell:
        r"""
        Rscript {input.script:q} {input.annotations:q} {output.source:q} {output.panel:q} \
            {output.png:q} {output.pdf:q} {input.summaries:q} > {log:q} 2>&1
        """
