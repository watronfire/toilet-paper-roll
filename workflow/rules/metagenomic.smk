rule kraken2_classify:
    input:
        reads1=lambda wildcards: SAMPLES[wildcards.sample]["read1"],
        reads2=lambda wildcards: SAMPLES[wildcards.sample]["read2"]
    params:
        db = KRAKEN_DB
    output:
        report = "intermediates/kraken2/{sample}_kraken2.txt"
    threads: 4
    shell:
        """
        kraken2 \
            --threads {threads} \
            -db {params.db} \
            --output - \
            --report {output.report} \
            --paired --gzip-compressed \
            {input.reads1} {input.reads2}
        """

rule normalize_taxonomy_reports:
    input:
        report = rules.kraken2_classify.output.report
    params:
        db = KRAKEN_DB,
        read_length = 150,
    output:
        corrected_report = "intermediates/kraken2/{sample}.bracken.txt",
        bracken_report = "intermediates/kraken2/{sample}.bracken.report"
    shell:
        """
        bracken \
            -d {params.db} \
            -r {params.read_length} \
            -i {input.report} \
            -o {output.bracken_report} \
            -w {output.corrected_report}
        """

rule combine_kraken2_reports:
    input:
        results = expand( "intermediates/kraken2/{sample}_kraken2.txt", sample=SAMPLES )
    output:
        summary = "results/reports/kraken2_reports.csv"
    run:
        toi = {
            "classified" : 1,
            "unclassified" : 0,
            "R2_bacteria" : 2,
            "R2_eurkaryote" : 2759,
            "R2_archaea" : 2157,
            "O_vibrionales" : 135623,
            "S_homo-sapiens" : 9606,
        }

        results = {
            "cell" : [],
        }
        for key in toi.keys():
            results[key] = []

        for f in input.results:
            sample = os.path.basename( f ).split( "_" )[0]
            tmp = pd.read_csv( f, sep="\t", header=None, names=["perc_fragments", "fragments_covered", "fragments_assigned", "rank", "ncbi_id", "sci_name" ] )
    
            results["cell"].append( sample )
            for k, v in toi.items():
                try:
                    results[k].append( tmp.loc[tmp["ncbi_id"]==v,"fragments_covered"].values[0] )
                except IndexError:
                    results[k].append( 0 )
        
        results = pd.DataFrame( results )
        results.to_csv( output.summary, index=False ) 


rule combine_bracken_reports:
    input:
        reports = expand( "intermediates/kraken2/{sample}.bracken.report", sample=SAMPLES )
    output:
        summary = "results/reports/bracken2_reports.csv"
    run:
        import pandas as pd
        
        bra = list()
        for f in input.reports:
            name = f.split("/")[-1].split( "." )[0]
            df = pd.read_csv( f, sep="\t" )
            df["cell"] = name
            bra.append( df )
        bra = pd.concat( bra )
        bra.to_csv( output.summary, index=False )
