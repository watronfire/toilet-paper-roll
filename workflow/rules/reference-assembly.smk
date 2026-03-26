rule alignment_bwa:
    message: "Mapping reads for {wildcards.sample} to {input.reference} using `bwa mem`."
    input:
        reads1=lambda wildcards: SAMPLES[wildcards.sample]["read1"],
        reads2=lambda wildcards: SAMPLES[wildcards.sample]["read2"],
        reference=REFERENCE,
        reference_index=rules.index_reference.output.reference_index
    params:
        bwa_params=config["alignment_bwa"]["bwa_params"]
    output:
        alignment="intermediates/illumina/merged_aligned_bams/{sample}.sorted.bam"
    threads: 8
    shell:
        """
        bwa mem \
            {params.bwa_params} \
            -t {threads} \
            {input.reference} \
            {input.reads1} {input.reads2} 2> /dev/null |\
        samtools view -Sb - |\
        samtools sort - |\
        samtools addreplacerg \
            -r "ID:{wildcards.sample}" \
            -o {output.alignment} - 
        """


# TODO: bedtools merge can be used to make depth_mask smaller.
rule generate_low_coverage_mask:
    message: "Create bed file from bam file for {wildcards.sample} indicating sites covered by less than {params.minimum_depth} reads"
    input:
        alignment=rules.alignment_bwa.output.alignment
    output:
        depth=temp( "intermediates/illumina/depth/{sample}.depth" ),
        depth_mask=temp( "intermediates/illumina/depth/{sample}.depthmask.bed" )
    params:
        minimum_depth=config["coverage_mask"]["required_depth"],
        minimum_base_quality=config["call_variants"]["minimum_base_quality"],
        minimum_mapping_quality=config["call_variants"]["minimum_mapping_quality"],
    shell:
        """
        samtools depth \
            -aa {input.alignment} \
            -q {params.minimum_base_quality} \
            -Q {params.minimum_mapping_quality} |\
        tee {output.depth} |\
        awk \
            -v depth="{params.minimum_depth}" \
            '$3 < depth {{printf "%s\\t%d\\t%d\\n", $1, $2 - 1, $2}}' \
            - > {output.depth_mask}
        """


rule call_variants_from_alignment:
    message: "Call variants from alignment for {wildcards.sample} using bcftools."
    input:
        alignment=rules.alignment_bwa.output.alignment,
        reference=REFERENCE,
        reference_index=rules.index_reference.output.reference_index
    params:
        maximum_depth=config["call_variants"]["maximum_depth"],
        minimum_mapping_quality=config["call_variants"]["minimum_mapping_quality"],
        minimum_base_quality=config["call_variants"]["minimum_base_quality"],
        mpileup_parameters=config["call_variants"]["mpileup_parameters"],
        call_parameters=config["call_variants"]["call_parameters"]
    output:
        variants="intermediates/illumina/variants/{sample}.bcftools.vcf"
    threads: 8
    shell:
        """
        bcftools mpileup \
            --threads {threads} \
            -d {params.maximum_depth} \
            -q {params.minimum_mapping_quality} \
            -Q {params.minimum_base_quality} \
            {params.mpileup_parameters} \
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
        filtered_variants="intermediates/illumina/variants/{sample}.bcftools.filt.vcf"
    group: "consensus"
    shell:
        """
        bcftools filter \
            --no-version \
            -i "INFO/AD[1]>{params.minimum_depth} && (INFO/AD[1])/(INFO/AD[0]+INFO/AD[1])>{params.minimum_support} && INFO/ADF[1]>{params.minimum_strand_depth} && INFO/ADR[1]>{params.minimum_strand_depth}" \
            -o {output.filtered_variants} \
            {input.variants}
        """

rule align_and_normalize_variants:
    message: "For sample {wildcards.sample}, Left-align and normalize indels, and remove insertions."
    input:
        variants=rules.filter_variants.output.filtered_variants,
        reference=REFERENCE,
        reference_index=rules.index_reference.output.reference_index
    output:
        normalized_variants="intermediates/illumina/variants/{sample}.bcftools.filt.norm.vcf.gz",
        variant_index="intermediates/illumina/variants/{sample}.bcftools.filt.norm.vcf.gz.csi"
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
        depth_mask=rules.generate_low_coverage_mask.output.depth_mask,
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

rule rename_fastq:
    message: "Generate temporary copy of FASTQ file for {wildcards.sample} to maintain names for QC."
    input:
        reads1=lambda wildcards: SAMPLES[wildcards.sample]["read1"]
    output:
        renamed_reads=temp( "intermediates/tmp/{sample}.fastq.gz" )
    shell:
        """
        cp {input.reads1} {output.renamed_reads}
        """


rule fastqc:
    # We're assuming that R1 is representative of R2. This should generally work and I can't think of a reason where
    # problems would only pop up in one rather than the other.
    message: "Calculate quality control metrics for raw sequencing reads of {wildcards.sample}."
    input:
        reads1=rules.rename_fastq.output.renamed_reads
    output:
        directory=directory( "results/reports/fastqc/{sample}/" ),
    threads: 2
    shell:
        """
        mkdir {output.directory} && \
        fastqc \
            --outdir {output.directory} \
            --threads {threads} \
            --quiet \
            {input.reads1}  
        """


rule alignment_stats:
    message: "Calculate the number of reads from {wildcards.sample} which map to the reference genome."
    input:
        alignment="intermediates/illumina/merged_aligned_bams/{sample}.sorted.bam"
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
        alignment="intermediates/illumina/merged_aligned_bams/{sample}.sorted.bam"
    output:
        reheaded_alignment="intermediates/illumina/merged_aligned_bams/{sample}.headed.bam",
        report_directory=directory( "results/reports/bamqc/{sample}/" )
    threads: 8
    shell:
        """
        samtools view -H {input.alignment} |\
        sed 's,^@RG.*,@RG\\tID:None\\tSM:None\\tLB:None\\tPL:Illumina,g' |\
        samtools reheader - {input.alignment} > {output.reheaded_alignment} && \
        qualimap bamqc \
            -bam {output.reheaded_alignment} \
            -nt {threads} \
            -outdir {output.report_directory}
        """



def get_qc_inputs( wildcards ):
    inputs = list()
    inputs.extend( expand( "results/reports/fastqc/{sample}/",sample=SAMPLES ) )
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

