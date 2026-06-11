rule fastq:
    message: "Calculate quality control metrics for raw sequencing reads of {wildcards.sample}."
    input:
        reads1=lambda wildcards: SAMPLES[wildcards.sample]["read1"],
        reads2=lambda wildcards: SAMPLES[wildcards.sample]["read2"]
    output:
        cleaned_r1 = "intermediates/trimmed_reads/{sample}_R1.trim.fastq.gz",
        cleaned_r2 = "intermediates/trimmed_reads/{sample}_R2.trim.fastq.gz",
        json_report="results/reports/fastqc/{sample}_fastp.json",
        html_report="results/reports/fastqc/{sample}_fastp.html",
    threads: 4
    group: "denovo"
    shell:
        """
        fastp \
            -i {input.reads1} -I {input.reads2} \
            -o {output.cleaned_r1} -O {output.cleaned_r2} \
            --json {output.json_report} \
            --html {output.html_report} \
            --report_title {wildcards.sample} \
            --thread {threads}
        """

rule megahit:
    input:
        r1 = rules.fastq.output.cleaned_r1,
        r2 = rules.fastq.output.cleaned_r2
    params:
        temp_outdir = "intermediates/assembly/{sample}/tmp"
    output:
        outdir = directory("intermediates/assembly/{sample}"),
        contigs = "intermediates/assembly/{sample}/{sample}.contigs.fa"
    threads: 8
    group: "denovo"
    shell:
        """
        rm -rf {params.temp_outdir}
        megahit -1 {input.r1} -2 {input.r2} -o {params.temp_outdir} --out-prefix {wildcards.sample} -t {threads} --memory 32000000000
        mv {params.temp_outdir}/* {output.outdir}
        """


rule classify_contigs:
    input:
        contigs = rules.megahit.output.contigs
    params:
        db = lambda wildcards: f"~/db/k2_standard-8"
    output:
        kraken_results = "results/kraken/{sample}.kraken.txt",
        kraken_report = "results/kraken/{sample}.report.txt"
    threads: 8
    group: "denovo"
    shell:
        """
        if [ ! -s "{input.contigs}" ]; then
            echo "Input contigs file is empty or missing. Creating empty placeholder outputs."
            touch "{output.kraken_results}" "{output.kraken_report}"
        else
            kraken2 --threads {threads} \
                --db {params.db} \
                --report {output.kraken_report} \
                --output {output.kraken_results} \
                {input.contigs}
        fi
        """


rule map_to_contigs:
    input:
        contigs = rules.megahit.output.contigs,
        reads1 = rules.fastq.output.cleaned_r1,
        reads2 = rules.fastq.output.cleaned_r2,
    output:
        alignment = "intermediates/assembly{sample}/{sample}.bam"
    threads: 8
    group: "denovo"
    shell:
        """
        minimap2 -ax sr \
            -t {threads} \
            {input.contigs} \
            {input.reads1} {input.reads2} |\
        samtools view -u -f2 - |\
        samtools sort -m 4G -o {output.alignment} -
        """


rule extract_vibrio_reads_from_contigs:
    input:
        alignment = rules.map_to_contigs.output.alignment,
        contigs = rules.megahit.output.contigs,
        kraken_results = rules.classify_contigs.output.kraken_results,
        kraken_report = rules.classify_contigs.output.kraken_report
    params:
        taxid = 135623,
        script = "workflow/scripts/extract_kraken_reads.py"
    output:
        vibrio_contigs = "intermediates/assembly/{sample}/{sample}.vibrio.fasta",
        vibrio_bed = "intermediates/assembly/{sample}/{sample}.vibrio.bed",
        filtered_reads1 = "intermediates/filtered_reads/{sample}.R1.fastq.gz",
        filtered_reads2 = "intermediates/filtered_reads/{sample}.R2.fastq.gz",
    group: "denovo"
    shell:
        """
        if [ ! -s "{input.kraken_results}" ]; then
            echo "Input results are empty or missing. Creating empty placeholder outputs."
            touch {output.vibrio_contigs} {output.vibrio_bed}
            echo -n "" | gzip > {output.filtered_reads1}
            echo -n "" | gzip > {output.filtered_reads2}
        else
            python {params.script} -k {input.kraken_results} -s {input.contigs} -t {params.taxid} -o {output.vibrio_contigs} --include-children -r {input.kraken_report} --noappend
        
        if [ -s "{output.vibrio_contigs}" ]; then
                fgrep ">" {output.vibrio_contigs} | sed "s/>//g" | awk '{{print $1"\t0\t1000000000"}}' > {output.vibrio_bed}
                samtools view -hb -L {output.vibrio_bed} {input.alignment} | samtools sort -m 4G -n - | samtools fastq -1 {output.filtered_reads1} -2 {output.filtered_reads2}
            else
                echo "No Vibrio reads found. Creating empty placeholder outputs."
                touch {output.vibrio_bed}
                echo -n "" | gzip > {output.filtered_reads1}
                echo -n "" | gzip > {output.filtered_reads2}
            fi
        fi
        """


rule alignment_minimap:
    message: "Mapping reads for {wildcards.sample} to {input.reference} using `bwa mem`."
    input:
        reads1=rules.extract_vibrio_reads_from_contigs.output.filtered_reads1,
        reads2=rules.extract_vibrio_reads_from_contigs.output.filtered_reads2,
        reference=REFERENCE,
        reference_index=rules.index_reference.output.reference_index
    output:
        alignment="intermediates/merged_aligned_bams/{sample}.sorted.bam"
    threads: 8
    group: "denovo"
    shell:
        """
        minimap2 -ax sr\
            -t {threads} \
            {input.reference} \
            {input.reads1} {input.reads2}  |\
        samtools view -hb -f2 - |\
        samtools sort -m 8G  - |\
        samtools addreplacerg \
            -r "ID:{wildcards.sample}" \
            -o {output.alignment} -
        """


# TODO: bedtools merge can be used to make depth_mask smaller.
rule generate_low_coverage_mask:
    message: "Create bed file from bam file for {wildcards.sample} indicating sites covered by less than {params.minimum_depth} reads"
    input:
        alignment=rules.alignment_minimap.output.alignment
    output:
        depth=temp( "intermediates/depth/{sample}.depth" ),
        depth_mask=temp( "intermediates/depth/{sample}.depthmask.bed" )
    params:
        minimum_depth=config["coverage_mask"]["required_depth"],
        minimum_base_quality=config["call_variants"]["minimum_base_quality"],
        minimum_mapping_quality=config["call_variants"]["minimum_mapping_quality"],
    group: "consensus"
    shell:
        """
        samtools view -u -f2 {input.alignment} |\
        samtools depth \
            -q {params.minimum_base_quality} \
            -Q {params.minimum_mapping_quality} \
            -aa - |\
        tee {output.depth} |\
        awk \
            -v depth="{params.minimum_depth}" \
            '$3 < depth {{printf "%s\\t%d\\t%d\\n", $1, $2 - 1, $2}}' \
            - > {output.depth_mask}
        """



rule call_variants_from_alignment:
    message: "Call variants from alignment for {wildcards.sample} using bcftools."
    input:
        alignment=rules.alignment_minimap.output.alignment,
        reference=REFERENCE,
        reference_index=rules.index_reference.output.reference_index
    params:
        maximum_depth=config["call_variants"]["maximum_depth"],
        minimum_mapping_quality=config["call_variants"]["minimum_mapping_quality"],
        minimum_base_quality=config["call_variants"]["minimum_base_quality"],
        mpileup_parameters=config["call_variants"]["mpileup_parameters"],
        call_parameters=config["call_variants"]["call_parameters"]
    output:
        variants="intermediates/variants/{sample}.bcftools.vcf"
    threads: 8
    shell:
        """
        bcftools mpileup \
            --threads {threads} \
            -d {params.maximum_depth} \
            -q {params.minimum_mapping_quality} \
            -Q {params.minimum_base_quality} \
            {params.mpileup_parameters} \
            --ignore-overlaps --skip-indels \
            -f {input.reference} \
            {input.alignment} |\
        bcftools call \
            --threads {threads} \
            {params.call_parameters} \
            -o {output.variants}
        """


rule filter_variants:
    message:
        """Remove variants for sample {wildcards.sample} that:
            - have depth less than {params.minimum_depth}
            - have individual strand depth less than {params.minimum_strand_depth}
            - are present in less than {params.minimum_support:.0%} of reads
        """
    input:
        variants=rules.call_variants_from_alignment.output.variants
    params:
        minimum_depth=config["filter_variants"]["minimum_depth"],
        minimum_strand_depth=config["filter_variants"]["minimum_strand_depth"],
        minimum_support=config["filter_variants"]["minimum_support"]
    output:
        filtered_variants="intermediates/variants/{sample}.bcftools.filt.vcf"
    group: "consensus"
    shell:
        """
        bcftools filter \
            --no-version \
            -i "INFO/AD[1]>={params.minimum_depth} && (INFO/AD[1])/(INFO/AD[0]+INFO/AD[1])>{params.minimum_support} && INFO/ADF[1]>={params.minimum_strand_depth} && INFO/ADR[1]>={params.minimum_strand_depth}" \
            -o {output.filtered_variants} \
            {input.variants}
        """


rule combine_depth_variants_mask:
    input:
        variants=rules.call_variants_from_alignment.output.variants,
        depth_mask=rules.generate_low_coverage_mask.output.depth_mask
    params:
        minimum_depth=config["filter_variants"]["minimum_depth"],
        minimum_strand_depth=config["filter_variants"]["minimum_strand_depth"],
        minimum_support=config["filter_variants"]["minimum_support"]
    output:
        bcftools_filtered_mask="intermediates/variants/{sample}.vcffiltered.bed",
        depth_filtered_mask="intermediates/variants/{sample}.allmask.bed"
    group: "consensus"
    shell:
        """
        bcftools filter \
            --no-version \
            -i "INFO/AD[1]<{params.minimum_depth} || (INFO/AD[1])/(INFO/AD[0]+INFO/AD[1])<={params.minimum_support} || INFO/ADF[1]<{params.minimum_strand_depth} || INFO/ADR[1]<{params.minimum_strand_depth}" \
            {input.variants} |\
        awk '(/^[^#]/ && length($4) == length($5)) {{printf "%s\\t%d\\t%d\\n", $1, $2 - 1, $2}}' > {output.bcftools_filtered_mask} &&\
        cat {output.bcftools_filtered_mask} {input.depth_mask} | sort -u > {output.depth_filtered_mask}
        """


rule align_and_normalize_variants:
    message: "For sample {wildcards.sample}, Left-align and normalize indels, and remove insertions."
    input:
        variants=rules.filter_variants.output.filtered_variants,
        reference=REFERENCE,
        reference_index=rules.index_reference.output.reference_index
    output:
        normalized_variants="intermediates/variants/{sample}.bcftools.filt.norm.vcf.gz",
        variant_index="intermediates/variants/{sample}.bcftools.filt.norm.vcf.gz.csi"
    group: "consensus"
    shell:
        """
        bcftools norm \
            --no-version \
            -f {input.reference} \
            {input.variants} |\
        bcftools filter \
            --no-version \
            --exclude 'strlen(REF)<strlen(ALT)' \
            -Oz \
            -o {output.normalized_variants} &&\
        bcftools index {output.normalized_variants}
        """


rule call_consensus:
    message: "For sample {wildcards.sample}, apply variants to reference to create consensus sequences. Masks sites with less than desired coverage."
    input:
        variants=rules.align_and_normalize_variants.output.normalized_variants,
        variant_index=rules.align_and_normalize_variants.output.variant_index,
        depth_mask=rules.combine_depth_variants_mask.output.depth_filtered_mask,
        reference=REFERENCE,
        reference_index=rules.index_reference.output.reference_index
    params:
        bcftools_parameters=config["call_consensus"]["consensus_parameters"]
    output:
        consensus_sequence="results/consensus/{sample}.consensus.fasta"
    group: "consensus"
    shell:
        """
        bcftools consensus \
            {params.bcftools_parameters} \
            --fasta-ref {input.reference} \
            --mask {input.depth_mask} \
            {input.variants} |\
        union -filter |\
        sed "1s/.*/>{wildcards.sample}/" > {output.consensus_sequence}
        """

def get_consensus_for_sample( wildcards ):
    return [f"results/consensus/{sample}.consensus.fasta" for sample in ORIGINAL_SAMPLES[wildcards.og]]


rule calculate_consensus_distance:
    input:
        consensus_sequences = lambda wildcards: expand( "results/consensus/{sample}.consensus.fasta", sample=ORIGINAL_SAMPLES[wildcards.og] )
    output:
        alignment = temp( "intermediates/tmp/{og}.alignment.fasta" ),
        distance = "intermediates/distances/{og}.distance.csv"
    shell:
        """
        cat {input.consensus_sequences} > {output.alignment} &&\
        pairsnp -sc {output.alignment} > {output.distance}
        """


rule combine_consensus_distances:
    input:
        distances = expand( "intermediates/distances/{og}.distance.csv", og=ORIGINAL_SAMPLES )
    output:
        distances = "results/reports/consensus_distances.csv"
    shell:
        """
        cat {input.distances} > {output.distances}
        """


rule alignment_stats:
    message: "Calculate the number of reads from {wildcards.sample} which map to the reference genome."
    input:
        alignment=rules.alignment_minimap.output.alignment
    output:
        alignment_stats="results/reports/samtools/{sample}.stats.txt",
        alignment_idxstats="results/reports/samtools/{sample}.idxstats.txt"
    shell:
        """
        samtools index {input.alignment} && \
        samtools idxstats {input.alignment} > {output.alignment_idxstats} && \
        samtools stats {input.alignment} > {output.alignment_stats}
        """


rule bamqc:
    message: "Assess the quality of the reference-based assembly of {wildcards.sample}."
    input:
        alignment=rules.alignment_minimap.output.alignment
    output:
        reheaded_alignment="intermediates/merged_aligned_bams/{sample}.headed.bam",
        report_directory=directory( "results/reports/bamqc/{sample}/" )
    threads: 8
    shell:
        """
        samtools reheader --command 'sed "s,^@RG.*,@RG\\tID:None\\tSM:None\\tLB:None\\tPL:Illumina,g"' {input.alignment} > {output.reheaded_alignment}
        if ! samtools view -F 0x4 -c {output.reheaded_alignment} | grep -q "^0$"; then
            qualimap bamqc \
                -bam {output.reheaded_alignment} \
                -nt {threads} \
                --java-mem-size=12G \
                -outdir {output.report_directory}
        else
            mkdir -p {output.report_directory}
            touch {output.report_directory}/genome_results.txt
        fi
        """


def get_qc_inputs( wildcards ):
    inputs = list()
    inputs.extend( expand( "results/reports/fastqc/{sample}_fastp.json",sample=SAMPLES ) )
    inputs.extend( expand( "results/reports/samtools/{sample}.stats.txt",sample=SAMPLES ) )
    inputs.extend( expand( "results/reports/samtools/{sample}.idxstats.txt",sample=SAMPLES ) )
    inputs.extend( expand( "results/reports/bamqc/{sample}/",sample=SAMPLES ) )
 #   inputs.extend( expand( "results/reports/quast/{sample}/", sample=SAMPLES ) )
    return inputs


rule generate_complete_report:
    message: "Combine individual QC reports into a single HTML report."
    input:
        get_qc_inputs
    params:
        multiqc_config=MULTIQC_CONFIG
    output:
        report="results/reports/qc_report.html",
        report_directory=directory( "results/reports/qc_report_data/" )
    shell:
        """
        multiqc \
            --filename {output.report} \
            --config {params.multiqc_config} \
            results/reports/
        """
