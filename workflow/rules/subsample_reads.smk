READS = [100_000, 250_000, 500_000, 1_000_000, 2_500_000, 5_000_000, "all"]
PBTREE = "~/scripts/vibecheck/vibecheck/resources/o1_cholera.no_missing.pb"


rule subsample_reads:
    input:
        read1 = lambda wildcards: SAMPLES[wildcards.sample]["read1"],
        bam = rules.alignment_minimap.output.alignment
    output:
        subsampled_names = "intermediates/subsamples/{sample}.{reads}.{trial}.txt",
        subsampled_bam = "intermediates/subsamples/{sample}.{reads}.{trial}.bam"
    group: "variants"
    run:
        if wildcards.reads != "all":
            reads = int( wildcards.reads )
            shell( """
                zcat {input.read1} | awk 'NR%4==1 {sub(/^@/, ""); print $1}' | shuf -n {wildcards.reads} > {output.subsampled_names} &&\
                samtools view -b -h -N {output.subsampled_names} {input.bam} > {output.subsampled_bam} &&\
                samtools index {output.subsampled_bam}
        """ )
        else:
            shell( "cp {input.bam} {output.subsampled_bam}" )
            shell( "touch {output.subsampled_names}" )


rule calculate_subsampled_depth:
    input:
        alignment = rules.subsample_reads.output.subsampled_bam
    params:
        minimum_depth=config["coverage_mask"]["required_depth"],
        minimum_base_quality=config["call_variants"]["minimum_base_quality"],
        minimum_mapping_quality=config["call_variants"]["minimum_mapping_quality"],
    output:
        depth = temp( "intermediates/subsampled_depth/{sample}.{reads}.{trial}.txt" ),
        coverage = "intermediates/subsampled_coverage/{sample}.{reads}.{trial}.txt"
    group: "variants"
    run:
        import pandas as pd

        shell( "samtools depth -aa -q {params.minimum_base_quality} -Q {params.minimum_mapping_quality} {input.alignment} > {output.depth}" )

        reads = wildcards.reads
        if wildcards.reads == "all":
            reads = shell( "samtools view -c {input.alignment}", read=True ).strip()

        df = pd.read_csv( output.depth, sep="\t", header=None, names=["ref", "pos", "depth"] )
        coverage = df.loc[df["depth"]>=params.minimum_depth].shape[0] / df.shape[0]
        depth = df["depth"].median()
        with open( output.coverage, "wt" ) as f:
            f.write( f"{wildcards.sample},{reads},{coverage},{depth}\n" )


rule get_variants_and_counts:
    input:
        alignment = rules.subsample_reads.output.subsampled_bam,
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
        variants = "intermediates/subsampled_variants/{sample}.{reads}.{trial}.variants.vcf.gz"
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
            --ignore-overlaps --skip-indels \
            -f {input.reference} \
            {input.alignment} |\
        tee >(bcftools query -f '[%CHROM\t%POS\t%REF\t%ALT\t%DP\t%AD]\n' - > {output.counts} ) |\
        bcftools call \
            --threads {threads} \
            {params.call_parameters} |\
        bcftools +fill-tags -Ou - -- -t AF |\
        bcftools filter \
            --no-version \
            -Oz \
            -i "INFO/AD[1]>={params.minimum_depth} && INFO/ADF[1]>={params.minimum_strand_depth} && INFO/ADR[1]>={params.minimum_strand_depth}" \
            -o {output.variants}
        bcftools index {output.variants}
        """


rule combine_vcfs:
    input:
        variants = expand( "intermediates/subsampled_variants/{sample}.{reads}.{trial}.variants.vcf.gz", sample=SAMPLES, reads=READS, trial=range(1) )
    output:
        combined_variants = "intermediates/usher/alignment.vcf"
    shell:
        """
        bcftools merge -Ov -o {output.combined_variants} {input.variants}
        """


rule cleanup_vcf:
    input:
        variants = rules.combine_vcfs.output.combined_variants
    output:
        cleaned_variants = "intermediates/usher/alignment.vcf.gz"
    run:
        import gzip
        import pandas as pd

        header_lines = []

        with open(input.variants, "rt") as f:
            for line in f:
                if line.startswith("#"):
                    header_lines.append(line)
                else:
                    # Once we hit data, we can stop reading line by line
                    break

        # The last header line is our column names (e.g., #CHROM, POS, ID...)
        # We need to strip the '#' from '#CHROM' so pandas can use it cleanly
        column_names = header_lines[-1].strip().split("\t")
        column_names[0] = "CHROM"

        # 2. Load the data into pandas, skipping the header rows
        # We use comment='#' to safely ignore the header lines we already captured
        df = pd.read_csv(
            input.variants,
            sep="\t",
            comment="#",
            header=None,
            names=column_names,
        )

        # 3. Modify the chromosome name and positions
        df.loc[df["CHROM"]=="AE003853","POS"] += 2961182
        df["CHROM"] = "AE003853"
        df = df.sort_values( "POS", ascending=True )

        # 4. Save the file back out as vcf.gz
        with gzip.open(output.cleaned_variants, "wt") as f_out:
            # Write the original header lines back first
            f_out.writelines(header_lines)

            # Write the pandas DataFrame data (excluding the pandas header)
            # We use lineterminator='\n' to ensure standard Unix line endings
            df.to_csv(f_out, sep="\t", index=False, header=False, lineterminator="\n")


rule classify_usher:
    input:
        variants = rules.cleanup_vcf.output.cleaned_variants
    params:
        pb_tree = PBTREE,
        outdir = "intermediates/vibecheck/"
    output:
        clades = "intermediates/vibecheck/clades.txt"
    threads: 8
    shell:
        """
        usher -n -D -i {params.pb_tree} -v {input.variants} -T {threads} -d {params.outdir}
        """


rule summarize_usher:
    input:
        clades = rules.classify_usher.output.clades
    output:
        results = "results/lineage_subsample.csv"
    run:
        import pandas as pd
        import re       
        import numpy as np

        vibe = { "sequence_file" : [], "lineage" : [], "confidence" : [], "classification_notes" : []}
        with open( input.clades, "rt" ) as f:
            for line in f:
                name, lineage_histogram = line.rstrip("\n").split("\t")
                if "*|" in lineage_histogram:
                    # example: A.28*|A.28(1/10),B.1(6/10),B.1.511(1/10),B.1.518(2/10)
                    lineage, histogram = lineage_histogram.split("*|")
                    histo_list = [i for i in histogram.split(",") if i]
                    confidence = 0.0
                    if len(histo_list) > 1:
                        for lin_counts in histo_list:
                            m = re.match(
                                r"([A-Z0-9.]+)\(([0-9]+)\/([0-9]+)\)", lin_counts
                            )
                            if m:
                                place_count, total = [
                                    int(m.group(2)),
                                    int(m.group(3)),
                                ]
                                confidence += (place_count / total) * np.log(
                                    place_count / total
                                )
                    confidence = np.exp(confidence)
                    histogram_note = "Usher placements: " + " ".join(histo_list)
                else:
                    lineage = lineage_histogram
                    confidence = 1.0
                    histogram_note = ""

                vibe["sequence_file"].append( name )
                vibe["lineage"].append( lineage )
                vibe["confidence"].append( confidence )
                vibe["classification_notes"].append( histogram_note)
        vibe = pd.DataFrame( vibe )
        vibe["cell"] = vibe["sequence_file"].str.extract( r"(plate\d-\w\d{1,2})" )
        vibe["reads"] = vibe["sequence_file"].str.extract( r".(\d+).0.bam")
        vibe = vibe.dropna( subset="reads")
        vibe["reads"] = vibe["reads"].astype(int)
        vibe["correct"] = ((vibe["lineage"]=="T12") & (vibe["confidence"]>0.5)).astype( int )
        vibe.to_csv( output.results, index=False )


rule summarize_variants_and_counts:
    input:
        alignment = rules.subsample_reads.output.subsampled_bam,
        counts = rules.get_variants_and_counts.output.counts,
        variants = rules.get_variants_and_counts.output.variants
    output:
        summary = "intermediates/subsampled_variants/{sample}.{reads}.{trial}.summary.txt"
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


        vcf = pd.read_csv( input.variants, sep="\t", header=None, names=["chrom", "pos", "id", "ref", "alt", "qual", "filter", "info", "format", "sample"],  comment="#", compression="gzip" )
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
        coverages = expand( "intermediates/subsampled_coverage/{sample}.{reads}.{trial}.txt", sample=SAMPLES, reads=READS, trial=range(1) )
    output:
        coverage_report = "results/coverage_subsample.csv"
    shell:
        """
        cat {input.coverages} > {output.coverage_report}
        """

rule combined_subsampled_variants:
    input:
        summaries = expand( "intermediates/subsampled_variants/{sample}.{reads}.{trial}.summary.txt", sample=SAMPLES, reads=[1_000_000], trial=range(1) )
    output:
        variant_report = "results/variants_subsample.csv"
    shell:
        """
        cat {input.summaries} > {output.variant_report}
        """
