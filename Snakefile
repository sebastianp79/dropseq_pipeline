#


import glob
import os

# Configuration ----------------------------------------------------------------

#expected number of cells (upper limit)
cell_num = config["cell_num"]

#cell barcode UMI structure
barcode = config["barcode"]

#genome_index
GenomeIndex = config["genome_index"]
#gene file
txn_file = config["txn_file"]

pd = config["proj_dir"]
data = pd + "data/"
output = pd + "output/"
fastq_dir = data + "fastq/"
fastqc_dir = output + "fastqc/"
fastq_merged = data + "fastq_merged/"
fastq_extr = data + "fastq_extr/"
fastq_trim = data + "fastq_trim/"
cell_stats = data + "cell_stats/"
aligned = data + "aligned/"
annotation_data = data + "annotation_data/"
sorted_reads = data + "sorted_reads/"
assigned = data + "assigned/"
dge_data = output + "dge_data/"
code = "Code/"


#make sure the project directory actually exists
assert os.path.exists(pd), "Project directory exists"

# Directory to send log files. Needs to be created manually since it
# is not a file created by a Snakemake rule.
dir_log = config["dir_log"]
if not os.path.isdir(dir_log):
    os.mkdir(dir_log)
# input data # might need to be changed to be universal

samples = set(glob_wildcards(fastq_dir + "{samples}_R1_001.fastq.gz").samples)


rule all:
    input:
        cell_stats + "whitelist.txt",
        fastq_extr + "combined_r2_extracted.fastq.gz",
        aligned + "Aligned.SortedByCoordinate.out.bam",
        sorted_reads + "assigned_sorted.bam",
        sorted_reads + "assigned_sorted.bam.bai",
        dge_data + "counts.tsv.gz",


rule run_fastqc:
    input:
        expand(fastqc_dir + "{sample}_R1_001_fastqc.html", sample = samples),
        expand(fastqc_dir + "{sample}_R2_001_fastqc.html", sample = samples)


rule fastqc:
    input:
        fastq_dir + "{sample}_R{read_num}_001.fastq.gz"
    output:
        fastqc_dir + "{sample}_R{read_num}_001_fastqc.html"
    params:
        outdir = fastqc_dir
    shell:
        "fastqc -o {params.outdir} {input}"


rule fastq_merge_r1:
    input:
        expand(fastq_dir + "{sample}_R1_001.fastq.gz", sample = samples)
    output:
        fastq_merged + "combined_r1.fastq.gz"
    shell:
        "zcat {input} | gzip --stdout > {output}"

rule fastq_merge_r2:
    input:
        expand(fastq_dir + "{sample}_R2_001.fastq.gz", sample = samples)
    output:
        fastq_merged + "combined_r2.fastq.gz"
    shell:
        "zcat {input} | gzip --stdout  > {output}"


# rule barcode_QC:
#     input:
#         r1 = fastq_merged + "combined_r1.fastq.gz",
#         r2 = fastq_merged + "combined_r2.fastq.gz"
#     output:
#         r1_trim = fastq_merged + "combined_r1_trimmed.fastq.gz",
#         r2_trim = fastq_merged + "combined_r2_trimmed.fastq.gz"
#     params:
#         min_len = 20,
#         min_qual = 20
#     shell:
#         "cutadapt --minimum-length {params.min_len} -q {params.min_qual}  -o {output.r1_trim} -p {output.r2_trim} {input.r1} {input.r1}"


rule umi_create_whitelist:
    input:
        fastq_merged + "combined_r1.fastq.gz"
    output:
        cell_stats + "whitelist.txt"
    params:
        cell_num = cell_num,
        bc = barcode
    shell:
        "umi_tools whitelist --stdin {input} --bc-pattern={params.bc} --plot-prefix=Whitelist_stats --extract-method=string --expect-cells={params.cell_num} --log2stderr > {output}"

rule umi_extract_bc_and_umi:
    input:
        r1 = fastq_merged + "combined_r1.fastq.gz",
        r2 = fastq_merged + "combined_r2.fastq.gz",
        wl = cell_stats + "whitelist.txt"
    output:
        r1_ext = temp(fastq_extr + "combined_r1_extracted.fastq.gz"),
        r2_ext = fastq_extr + "combined_r2_extracted.fastq.gz"
    params:
        bc = barcode
    shell:
        "umi_tools extract --bc-pattern={params.bc} --stdin {input.r1} --stdout {output.r1_ext} --read2-in {input.r2} --read2-out={output.r2_ext}  --error-correct-cell --filter-cell-barcode --whitelist={input.wl}"


#allows for proportion of 0.1 mismatches (e.g. 10*0.1=1 mismatch in 10bp adaptor match, bp below zero are rounded to 0)
# sequences to match are poly at the 3 prime end (at least 6) and TSO oligo on the 5 prime end
rule trim_read2:
    input:
        r1 = fastq_extr + "combined_r2_extracted.fastq.gz"
    output:
        r1_trim = fastq_extr + "combined_r2_extracted_polyA_adaptor_trimmed.fastq.gz",
    shell:
        "cutadapt -a AAAAAA -g AAGCAGTGGTATCAACGCAGAGTGAATGGG -o {output} {input}"


rule align:
    input:
        fq = fastq_extr + "combined_r2_extracted_polyA_adaptor_trimmed.fastq.gz",
        ref_genome = GenomeIndex
    output:
        aligned + "Aligned.SortedByCoordinate.out.bam"
    threads: 4
    shell:
        "STAR --runThreadN {threads} --genomeDir {input.ref_genome} --readFilesIn {input.fq} --readFilesCommand zcat --outFilterMultimapNmax 1 --outSAMtype BAM SortedByCoordinate --outStd BAM_SortedByCoordinate > {output}"

# assign reads to transcripts

rule reads_to_transcripts:
    input:
        bam = aligned + "Aligned.SortedByCoordinate.out.bam",
        features = txn_file
    output:
        assigned_feat = assigned + "gene_assigned.summary",
        bam_counts = assigned + "Aligned.SortedByCoordinate.out.bam.featureCounts.bam"
    threads: 4
    shell:
        "featureCounts -a {input.features} -o {output.assigned_feat} -R BAM {input.bam} -T {threads}"

#might need to change -m memory allocation per thread
rule sort_bams:
    input:
        assigned + "Aligned.SortedByCoordinate.out.bam.featureCounts.bam"
    output:
        sorted_reads + "assigned_sorted.bam"
    shell:
        "samtools sort -o {output} -O bam {input}"

rule index_bams:
    input:
        sorted_reads + "assigned_sorted.bam"
    output:
        sorted_reads + "assigned_sorted.bam.bai"
    shell:
        "samtools index {input}"

rule make_DGE_matrix:
    input:
        sorted_reads + "assigned_sorted.bam",
        sorted_reads + "assigned_sorted.bam.bai"
    output:
        dge_data + "counts.tsv.gz"
    shell:
        "umi_tools count --wide-format-cell-counts --per-gene --gene-tag=XT --per-cell -I {input} -S {output}"

# rule qc_with_multiqc:
#     input:
#         fqc1 = not sure how
#         fqc2 = alternative input dir,
#         alternative input dir
#         alternative input dir
#     output:
#         multiqc_report.html
#     shell:
#         "multiqc "
#         multiqc data/
# multiqc data/ ../proj_one/analysis/ /tmp/results
# multiqc data/*_fastqc.zip
# multiqc data/sample_1*
#
#
#
#
# The report is called multiqc_report.html by default. Tab-delimited data files are created in multiqc_data/, containing additional information. You can use a custom name for the report with the -n/--filename parameter, or instruct MultiQC to create them in a subdirectory using the -o/-outdir parameter.
#
# Note that different MultiQC templates may have different defaults.
