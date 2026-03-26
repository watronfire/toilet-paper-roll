READS = [1_000, 10_000, 100_000, 1_000_000, 5_000_000, "all"]

rule subsample_reads:
    input:
        alignment = rules.alignment_bwa.output.alignment
    output:
        subalignment = temp( "intermediates/subsamples/{sample}.{reads}.{trial}.bam" )
    run:
        if wildcards.reads != "all":
            reads = int( wildcards.reads )        
            shell( """
        cat <(samtools view -H {input.alignment}) <(samtools view {input.alignment} | shuf -n {reads}) |\
        samtools view -b - |\
        samtools sort - > {output.subalignment}
        """ )
        else:
            shell( "cp {input.alignment} {output.subalignment}" )


rule calculate_subsampled_depth:
    input:
        alignment = rules.subsample_reads.output.subalignment
    params:
        minimum_depth=config["coverage_mask"]["required_depth"],
        minimum_base_quality=config["call_variants"]["minimum_base_quality"],
        minimum_mapping_quality=config["call_variants"]["minimum_mapping_quality"],
    output:
        depth = temp( "intermediates/subsamples/{sample}.{reads}.{trials}.txt" ),
        coverage = "intermediates/subsampled_coverage/{sample}.{reads}.{trials}.txt"
    run:
        import pandas as pd

        shell( "samtools depth -aa -q {params.minimum_base_quality} -Q {params.minimum_mapping_quality} {input.alignment} > {output.depth}" )
        
        reads = wildcards.reads
        if wildcards.reads == "all":
            reads = shell( "samtools view -c {input.alignment}", read=True ).strip()
        
        df = pd.read_csv( output.depth, sep="\t", header=None, names=["ref", "pos", "depth"] )
        coverage = df.loc[df["depth"]>params.minimum_depth].shape[0] / df.shape[0]
        depth = df["depth"].median()
        with open( output.coverage, "wt" ) as f:
            f.write( f"{wildcards.sample},{reads},{coverage},{depth}\n" )


rule get_variants_and_counts:
    input:
        alignment = rules.subsample_reads.output.subalignment,
        reference = REFERENCE
    params:
        maximum_depth=config["call_variants"]["maximum_depth"],
        minimum_mapping_quality=config["call_variants"]["minimum_mapping_quality"],
        minimum_base_quality=config["call_variants"]["minimum_base_quality"],
        mpileup_parameters=config["call_variants"]["mpileup_parameters"],
        call_parameters="-mv -Ou --ploidy 1",
        minimum_depth=config["filter_variants"]["minimum_depth"],
        minimum_strand_depth=config["filter_variants"]["minimum_strand_depth"],
    output:
        counts = "intermediates/subsampled_variants/{sample}.{reads}.{trial}.counts.txt",
        variants = "intermediates/subsampled_variants/{sample}.{reads}.{trial}.variants.vcf"
    threads: 4
    group: "variants"
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
        tee >(bcftools query -f '[%CHROM\t%POS\t%REF\t%ALT\t%DP\t%AD]\n' - > {output.counts} ) |\
        bcftools call \
            --threads {threads} \
            {params.call_parameters} |\
        bcftools +fill-tags -Ou - -- -t AF |\
        bcftools filter \
            --no-version \
            -Ov \
            -i "INFO/AD[1]>{params.minimum_depth} && INFO/ADF[1]>{params.minimum_strand_depth} && INFO/ADR[1]>{params.minimum_strand_depth}" \
            -o {output.variants}
        """

rule summarize_variants_and_counts:
    input:
        alignment = rules.subsample_reads.output.subalignment,
        counts = rules.get_variants_and_counts.output.counts,
        variants = rules.get_variants_and_counts.output.variants
    output:
        summary = "intermediates/subsample_variants/{sample}.{reads}.{trial}.summary.txt"
    group: "variants"
    run:
        import pandas as pd
        import numpy as np

        reads = shell( "samtools view -c {input.alignment}", read=True ).strip()

        counts = pd.read_csv( input.counts, sep="\t", header=None, names=["chrom", "pos", "ref", "alt", "depth", "counts"] )
        if counts.shape[0] == 0:
            pi = None
            pi_median = None
        else:
            counts["diversity"] = counts.apply( lambda x:  1 - np.sum( np.power(np.array( list(map( int, x["counts"].split( "," ) ) ) ) / x["depth"], 2) ), axis=1 )
            pi = counts["diversity"].mean()
            pi_median = counts["diversity"].median()


        vcf = pd.read_csv( input.variants, sep="\t", header=None, names=["chrom", "pos", "id", "ref", "alt", "qual", "filter", "info", "format", "sample"],  comment="#" )
        vcf = vcf.loc[vcf["ref"].str.len()==vcf["alt"].str.len()]
        vcf["AD"] = vcf["info"].str.extract( r"AD=([0-9,]+);" )
        vcf["AD"] = vcf["AD"].str.split( "," )
        assert (vcf["AD"].str.len() == 2).all()
        vcf["AF"] = vcf["AD"].apply( lambda x: int( x[1] ) / (sum( map( int, x ) ) ) )
        vcf = vcf.loc[(vcf["AF"] > 0.03)&(vcf["AF"]<1)]
        vcf.loc[vcf["AF"] > 0.5, "AF"] = 1 - vcf["AF"]
        af = vcf["AF"].mean()
        af_median = vcf["AF"].median()
        variant_count = vcf.shape[0]

        with open( output.summary, "wt" ) as outf:
            outf.write( f"{wildcards.sample},{reads},{pi},{pi_median},{af},{af_median},{variant_count}\n" )


rule combine_subsampled_coverage:
    input:
        coverages = expand( "intermediates/subsampled_coverage/{sample}.{reads}.{trials}.txt", sample=SAMPLES, reads=READS, trials=range(3) )
    output:
        coverage_report = "results/coverage_subsample.csv" 
    shell:
        """
        cat {input.coverages} > {output.coverage_report}
        """

rule combined_subsampled_variants:
    input:
        summaries = expand( "intermediates/subsample_variants/{sample}.{reads}.{trial}.summary.txt", sample=SAMPLES, reads=[1_000_000], trial=range(1) )
    output:
        variant_report = "results/variants_subsample.csv"
    shell:
        """
        cat {input.summaries} > {output.variant_report}
        """

