rule quality_trimming:
    message:
        """Trim low quality bases from raw sequencing reads of {wildcards.sample}, using fastp.
        """
    input:
        reads1=lambda wildcards: SAMPLES[wildcards.sample]["read1"],
        reads2=lambda wildcards: SAMPLES[wildcards.sample]["read2"]
    output:
        read1_trimmed="intermediates/illumina/trimmed/{sample}_R1.fastq.gz",
        read2_trimmed="intermediates/illumina/trimmed/{sample}_R2.fastq.gz",
        unpaired=temp( "intermediates/illumina/trimmed/{sample}_U.fastq.gz" ),
        json_report=temp( "intermediates/illumina/trimmed/{sample}.json" ),
        html_report=temp( "intermediates/illumina/trimmed/{sample}.html" )
    group: "denovo_{sample}"
    threads: 8
    shell:
        """
        fastp \
            --in1 {input.reads1} \
            --in2 {input.reads2} \
            --out1 {output.read1_trimmed} \
            --out2 {output.read2_trimmed} \
            --unpaired1 {output.unpaired} \
            --unpaired2 {output.unpaired} \
            -h {output.html_report} \
            -j {output.json_report}
        """


rule estimate_genome_size:
    input:
        reads=rules.quality_trimming.output.read1_trimmed
    params:
        prefix="intermediates/illumina/{sample}/foo",
        kmer_length=21,
        minimum_occurance=10
    output:
        estimates="intermediates/illumina/trimmed/{sample}.txt",
        intermediates=temp( directory( "intermediates/illumina/temporary/{sample}/" ) )
    threads: 8
    resources:
        mem_gb=8
    group: "denovo_{sample}"
    shell:
        """
        mkdir {output.intermediates} &&\
        kmc \
            -sm -w \
            -m{resources.mem_gb} \
            -t{threads} \
            -k{params.kmer_length} \
            -ci{params.minimum_occurance} \
            {input.reads} \
            {params.prefix} \
            {output.intermediates} > {output.estimates}
        """


rule read_correction:
    message: "Corrects sequencing errors in the reads for {wildcards.sample} by sampling kmers with lighter."
    input:
        reads1=rules.quality_trimming.output.read1_trimmed,
        reads2=rules.quality_trimming.output.read2_trimmed,
        genome_stats=rules.estimate_genome_size.output.estimates
    params:
        kmer_length=32,
        maximum_corrections=1,
        output_directory="intermediates/illumina/corrected_reads/"
    output:
        corrected_reads1="intermediates/illumina/corrected_reads/{sample}_R1.cor.fq.gz",
        corrected_reads2="intermediates/illumina/corrected_reads/{sample}_R2.cor.fq.gz"
    threads: 8
    group: "denovo_{sample}"
    shell:
        """
        GENOMESIZE=$(fgrep "unique counted" {input.genome_stats} | grep -Eo "(\S+)\s*$" | tail -n1) &&\
        lighter \
            -r {input.reads1} \
            -r {input.reads2} \
            -K {params.kmer_length} $GENOMESIZE \
            -t {threads} \
            -maxcor {params.maximum_corrections} \
            -od {params.output_directory}
        """


rule read_stitching:
    message: "Stitch overlapping reads together in {wildcards.sample} using Flash."
    input:
        reads1=rules.read_correction.output.corrected_reads1,
        reads2=rules.read_correction.output.corrected_reads2,
        raw_reads1=lambda wildcards: SAMPLES[wildcards.sample]["read1"]
    params:
        output_directory="intermediates/illumina/stitched_reads/",
        minimum_overlap=20
    output:
        stitched_reads="intermediates/illumina/stitched_reads/{sample}.extendedFrags.fastq.gz",
        unstitched_read1="intermediates/illumina/stitched_reads/{sample}.notCombined_1.fastq.gz",
        unstitched_read2="intermediates/illumina/stitched_reads/{sample}.notCombined_2.fastq.gz",
        histogram=temp( "intermediates/illumina/stitched_reads/{sample}.hist" ),
        histogram2=temp( "intermediates/illumina/stitched_reads/{sample}.histogram" )
    threads: 8
    group: "denovo_{sample}"
    shell:
        """
        MAXOVERLAP=$(zcat < {input.raw_reads1} | awk '{{if(NR%4==2) {{count++; bases += length}} }} END{{print int( bases/count )}}') &&\
        flash \
            -m {params.minimum_overlap} \
            -M $MAXOVERLAP \
            -d {params.output_directory} \
            -o {wildcards.sample} \
            -z \
            -t {threads} \
            {input.reads1} \
            {input.reads2}
        """


rule denovo_assembly:
    message: "Assemble reads using a wrapper around SPAdes for {wildcards.sample}"
    input:
        reads1=rules.read_stitching.output.unstitched_read1,
        reads2=rules.read_stitching.output.unstitched_read2,
        stitched_reads=rules.read_stitching.output.stitched_reads
    params:
        unicycle_params="--min_component_size 200 --keep 0",
        temp_assembly="intermediates/illumina/assembly/{sample}/assembly.fasta",
        temp_graph="intermediates/illumina/assembly/{sample}/assembly.gfa",
        temp_log="intermediates/illumina/assembly/{sample}/unicycler.log"
    output:
        assembly="intermediates/illumina/assembly/{sample}.assembly.fasta",
        temporary_directory=temp( directory( "intermediates/illumina/assembly/{sample}/" ) ),
        assembly_graph="intermediates/illumina/assembly/{sample}.gfa",
        assembly_log="intermediates/illumina/assembly/{sample}.log"
    threads: 8
    group: "denovo_{sample}"
    shell:
        """
        unicycler \
            -1 {input.reads1} \
            -2 {input.reads2} \
            -s {input.stitched_reads} \
            -o {output.temporary_directory} \
            --threads {threads} \
            {params.unicycle_params} &&\
        mv {params.temp_assembly} {output.assembly} &&\
        mv {params.temp_graph} {output.assembly_graph} &&\
        mv {params.temp_log} {output.assembly_log}
        """


rule assembly_stats:
    message: "Assess the quality of the de novo assembly for {wildcards.sample}."
    input:
        assembly=rules.denovo_assembly.output.assembly
    params:
        quast_arguments="--fast --space-efficient"
    output:
        report_directory=directory( "results/reports/quast/{sample}/" )
    threads: 8
    group: "denovo_{sample}"
    shell:
        """
        quast \
            {params.quast_arguments} \
            -o {output.report_directory} \
            {input.assembly}
        """
