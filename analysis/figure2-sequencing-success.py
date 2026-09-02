import marimo

__generated_with = "0.23.16"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo
    import pandas as pd
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    import matplotlib.ticker as mticker
    from utils import setup_plotting_standards, basic_formatting, get_sample_types
    import numpy as np

    from statsmodels.stats.proportion import proportion_confint

    setup_plotting_standards()

    sample_type_order, sample_type_long_dict, sample_type_display_dict = get_sample_types()

    samples = pd.read_csv( "resources/samples.csv" )
    return (
        basic_formatting,
        mticker,
        np,
        pd,
        plt,
        proportion_confint,
        sample_type_display_dict,
        sample_type_order,
        samples,
    )


@app.cell
def _(pd, proportion_confint, samples):
    rare = pd.read_csv( "results/coverage_subsample.csv", header=None, names=["cell","reads", "coverage", "median_depth"])
    rare = rare.merge( samples, on="cell", how="left" )
    success_per_reads = rare.groupby( ["reads","sample_type"] ).agg( success=("coverage", lambda x: (x > .9).sum()), count=("coverage", "count") )
    success_per_reads["prob"] = success_per_reads["success"] / success_per_reads["count"]
    success_per_reads["lower"], success_per_reads["upper"] = proportion_confint(count=success_per_reads["success"], nobs=success_per_reads["count"], alpha=0.05, method='jeffreys')
    success_per_reads["err_low"] = (success_per_reads["prob"] - success_per_reads["lower"]).abs()
    success_per_reads["err_high"] = (success_per_reads["upper"] - success_per_reads["prob"]).abs()
    success_per_reads
    return (success_per_reads,)


@app.cell
def _(
    basic_formatting,
    mticker,
    np,
    pd,
    plt,
    sample_type_display_dict,
    sample_type_order,
    success_per_reads,
):
    cmap = {"Isolates" : "#66c2a5", "APW" : "#fc8d62", "Stool" : "#8da0cb"}

    bar_props = { 
        "zorder" : 10, 
        "width" : 0.85, 
        "edgecolor" : "black", 
        "linewidth" : 0.75
    }
    error_props = {
        "zorder" : 15, 
        "fmt" : "none", 
        "capsize" : 3, 
        "capthick" : 0.75, 
        "elinewidth" : 0.75, 
        "ecolor" : "black"
    }

    col_config = [
        (0, sample_type_order[:4], "Isolates"),
        (1, sample_type_order[4:6], "APW"),
        (2, sample_type_order[6:], "Stool"),
    ]

    for _num in [1_000_000, 5_000_000]:
        _fig, _axes = plt.subplots( figsize=(6,4), ncols=3, gridspec_kw={"width_ratios" : [4, 2, 4]})

        for _ax, (_col_idx, _sample_subset, _x_label) in zip( _axes, col_config ):
            plot_df = success_per_reads.loc[pd.IndexSlice[_num, _sample_subset], :]
            _ax.bar( x=range( len( _sample_subset ) ), height=plot_df["prob"], color=cmap[_x_label], **bar_props )
            _ax.errorbar( range( len( _sample_subset ) ), plot_df["prob"], yerr=plot_df[["err_low", "err_high"]].T, **error_props )

            _ax.set_xticks( range( len( _sample_subset ) ), [sample_type_display_dict[i].split( "\n" )[0] for i in _sample_subset] )

            _ax.set_yticks(np.arange(0,1.2,0.2))
            _ax.set_yticks(np.arange(0,1.0,0.05), minor=True )
            _ax.yaxis.set_major_formatter( mticker.PercentFormatter(1) )

            basic_formatting( _ax, which="y", ylabel="Proportion with complete genome", ylims=(-0.01,1.02), xsize=10, ysize=12, spines=[] )
            _ax.set_xlabel( _x_label, fontsize=10, fontweight="medium")
            if _col_idx == 0:
                _ax.set_title( f"{_num:,} reads", fontsize=12, fontweight="medium", loc="left" )

            _ax.label_outer()

        plt.tight_layout(w_pad=0.5)
        plt.savefig( f"analysis/plots/figure2-sequencing-success-{_num:.0f}.pdf" )
        plt.show()
    return bar_props, cmap, col_config, error_props


@app.cell
def _(pd, proportion_confint, samples):
    line = pd.read_csv( "results/lineage_subsample.csv" )
    line = line.merge( samples, on="cell", how="left" )
    line = line.loc[line["reads"].isin([100_000,250_000,500_000,1_000_000,2_500_000,5_000_000])].dropna( subset=["correct"] )

    line_per_reads = line.groupby( ["reads","sample_type"] ).agg( success=("correct", "sum"), count=("correct", "count") )
    line_per_reads["prob"] = line_per_reads["success"] / line_per_reads["count"]
    line_per_reads["lower"], line_per_reads["upper"] = proportion_confint(count=line_per_reads["success"], nobs=line_per_reads["count"], alpha=0.05, method='jeffreys')
    line_per_reads["err_low"] = (line_per_reads["prob"] - line_per_reads["lower"]).abs()
    line_per_reads["err_high"] = (line_per_reads["upper"] - line_per_reads["prob"]).abs()
    line_per_reads
    return (line_per_reads,)


@app.cell
def _(
    bar_props,
    basic_formatting,
    cmap,
    col_config,
    error_props,
    line_per_reads,
    mticker,
    np,
    pd,
    plt,
    sample_type_display_dict,
):
    for _num in [1_000_000, 5_000_000]:
        _fig, _axes = plt.subplots( figsize=(6,4), ncols=3, gridspec_kw={"width_ratios" : [4, 2, 4]})

        for _ax, (_col_idx, _sample_subset, _x_label) in zip( _axes, col_config ):
            _plot_df = line_per_reads.loc[pd.IndexSlice[_num, _sample_subset], :]
            _ax.bar( x=range( len( _sample_subset ) ), height=_plot_df["prob"], color=cmap[_x_label], **bar_props )
            _ax.errorbar( range( len( _sample_subset ) ), _plot_df["prob"], yerr=_plot_df[["err_low", "err_high"]].T, **error_props )

            _ax.set_xticks( range( len( _sample_subset ) ), [sample_type_display_dict[i].split( "\n" )[0] for i in _sample_subset] )

            _ax.set_yticks(np.arange(0,1.2,0.2))
            _ax.set_yticks(np.arange(0,1.0,0.05), minor=True )
            _ax.yaxis.set_major_formatter( mticker.PercentFormatter(1) )

            basic_formatting( _ax, which="y", ylabel="Proportion with accurate lineage", ylims=(-0.01,1.02), xsize=10, ysize=12, spines=[] )
            _ax.set_xlabel( _x_label, fontsize=10, fontweight="medium")
            if _col_idx == 0:
                _ax.set_title( f"{_num:,} reads", fontsize=12, fontweight="medium", loc="left" )

            _ax.label_outer()

        plt.tight_layout(w_pad=0.5)
        plt.savefig( f"analysis/plots/figure2-classification-success-{_num:.0f}.pdf" )
        plt.show()
    return


if __name__ == "__main__":
    app.run()
