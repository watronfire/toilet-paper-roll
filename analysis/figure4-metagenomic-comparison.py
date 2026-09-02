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
    import seaborn as sns
    from scipy.spatial import distance

    setup_plotting_standards()

    sample_type_order, sample_type_long_dict, sample_type_display_dict = get_sample_types()

    samples = pd.read_csv( "resources/samples.csv" )
    return (
        Line2D,
        basic_formatting,
        distance,
        np,
        pd,
        plt,
        sample_type_display_dict,
        sample_type_order,
        samples,
        sns,
    )


@app.cell
def _(pd, samples):
    bra = pd.read_csv( "results/reports/bracken2_reports.csv" )
    bra = bra.merge( samples, on="cell", how="left" )
    return (bra,)


@app.cell
def _(bra, np, samples):
    shannon = bra.groupby( "cell" ).apply( lambda x: -np.sum(x["fraction_total_reads"] *  np.log( x["fraction_total_reads"])) )
    shannon = shannon.reset_index().rename( columns={0:"shannon"} )
    shannon = shannon.merge( samples, on="cell" )
    shannon
    return (shannon,)


@app.cell
def _(
    basic_formatting,
    np,
    plt,
    sample_type_display_dict,
    sample_type_order,
    shannon,
    sns,
):
    # 1. Store common boxplot styling in a single dictionary
    _cmap = {"Isolates" : "#66c2a5", "APW" : "#fc8d62", "Stool" : "#8da0cb"}
    _box_props = { "markeredgecolor": "black", "markeredgewidth": 0.75 }
    _style_kwargs = {
        "linecolor": "black",
        "fliersize": 5,
        "linewidth": 0.75,
        "formatter": lambda x: sample_type_display_dict[x].split("\n")[0],
        "medianprops": {"linewidth": 1.5},
        "capwidths": 0.15,
    }

    # 3. Create subplots
    _fig, _axes = plt.subplots(
        figsize=(8, 4), ncols=3, gridspec_kw={"width_ratios": [4, 2, 4]}, sharey=True
    )

    # 2. Define panel grouping configurations
    panel_config = [
        (sample_type_order[:4], "Isolates"),
        (sample_type_order[4:6], "APW"),
        (sample_type_order[6:], "Stool"),
    ]


    # 4. Loop through panels to plot data and apply formatting
    for _ax, (_sample_subset, _label) in zip( _axes, panel_config ):
        _subset_data = shannon.loc[shannon["sample_type"].isin(_sample_subset)]

        sns.boxplot(
            data=_subset_data,
            x="sample_type",
            y="shannon",
            order=_sample_subset,
            ax=_ax,
            color=_cmap[_label],
            flierprops={ "markeredgecolor": "black", "markeredgewidth": 0.75, "markerfacecolor" : _cmap[_label]},
            **_style_kwargs,
        )

        _ax.set_yticks( np.arange(0,6,1), )
        _ax.set_yticks( np.arange(0,5,0.25), minor=True )

        basic_formatting(
            _ax,
            ylabel="Shannon index\n(Higher is more diverse)",
            ylims=(-0.1, 5.1),
            xsize=14,
            ysize=14,
            spines=[],
        )
    
        _ax.tick_params(axis="both", which="both", size=0)
        _ax.label_outer()
        _ax.set_xlabel( _label, fontsize=14, fontweight="medium" )

    plt.tight_layout(w_pad=0.5)
    plt.savefig( "analysis/plots/figure4-shannon-diversity.pdf" )
    plt.show()
    return


@app.cell
def _(basic_formatting, bra, distance, np, pd, plt, sns):
    _box_props = { "markeredgecolor": "black", "markeredgewidth": 0.75, "markerfacecolor" : "#E5E4E2" }
    _style_kwargs = {
        "linecolor": "black",
        "fliersize": 5,
        "linewidth": 0.75,
        "formatter": lambda x: x,
        "medianprops": {"linewidth": 1.5},
        "capwidths": 0.15,
        "color" : "#E5E4E2"
    }

    fig, ax= plt.subplots( figsize=(8,4) )
    all_samples = bra["sample"].sort_values().unique()
    all_bc = list()
    for idx, (sample, df) in enumerate( bra.groupby( "sample" ) ):
        sample_df = bra.loc[bra["sample"]==sample].pivot_table( index="name", columns="sample_type", values="new_est_reads", fill_value=0 )
        sample_df = sample_df.drop(index=["Homo sapiens", "Vibrio cholerae"] )
        bc = distance.pdist( sample_df.T, metric="braycurtis" )

        points= sns.boxplot( x=idx, y=bc, ax=ax, zorder=10, flierprops=_box_props, **_style_kwargs )
        points.set_clip_on(False)

        bc = pd.DataFrame( distance.squareform( bc ), index=sample_df.columns, columns=sample_df.columns ).melt(ignore_index=False, var_name="stB", value_name=sample ).reset_index().rename( columns={"sample_type" : "stA"}).set_index( ["stA", "stB"] )
        all_bc.append( bc )

    all_bc = pd.concat( all_bc, ignore_index=False, axis=1 )

    ax.set_xticks(range(0, len(all_samples)), [int( i[2:] ) for i in all_samples] )
    ax.set_yticks( np.arange(0,1.2, 0.2) )
    ax.set_yticks( np.arange(0,1, 0.05), minor=True )

    basic_formatting(
        ax, 
        xlabel="Sample", 
        ylabel="Bray-Curtis distance\n(Lower is more similar)", 
        xlims=(-0.5,len(all_samples)-0.5), 
        ylims=(-0.01,1.01), 
        spines=[],
        xsize=14,
        ysize=14,
    )

    plt.tight_layout()
    plt.savefig( "analysis/plots/figure4-distance.pdf" )
    plt.show()
    return all_bc, all_samples


@app.cell
def _(all_bc, plt, sample_type_order, sns):
    sns.heatmap( all_bc.mean( axis=1 ).reset_index().rename(columns={0 : "mean_dist" } ).pivot_table( index="stA", columns="stB", values="mean_dist" ).reindex( index=sample_type_order, columns=sample_type_order ), vmin=0, vmax=1, cmap="magma_r" )
    #plt.tight_layout()
    plt.savefig( "analysis/plots/unknown-bray-curtis-distance.pdf" )
    plt.show()
    return


@app.cell
def _(plt):
    plt.get_cmap()
    return


@app.cell
def _(
    Line2D,
    all_samples,
    bra,
    pd,
    plt,
    sample_type_display_dict,
    sample_type_order,
    samples,
):
    cmap = list( plt.get_cmap("tab20").colors )
    cmap.append( (0.86,0.86,0.86) )

    pt_reads = bra.pivot_table( index="name", columns="cell", values="new_est_reads", fill_value=0 )
    pt_reads = pt_reads.drop(index=["Homo sapiens", "Vibrio cholerae"] )

    top20 = (pt_reads / pt_reads.sum()).mean(axis=1).sort_values().tail(20).index

    other = pt_reads.loc[[i for i in pt_reads.index if i not in top20]].sum()
    other.name = "Other"

    pt_reads = pd.concat( [pt_reads.T[top20], other.T], axis=1 )

    _fig, _axes = plt.subplots( figsize=(14, 6.5), nrows=len(sample_type_order), ncols=len( all_samples ) )

    for row, st in zip( _axes, sample_type_order ):
        for _ax, _sample in zip( row, all_samples):
            try:
                cell = samples.loc[(samples["sample"]==_sample)&(samples["sample_type"]==st),"cell"].values[0]
            except IndexError:
                [_ax.spines[j].set_visible( False ) for j in _ax.spines]
                _ax.set_yticks([])
                _ax.set_xticks([])
                _ax.set_ylim(-1,1)
                _ax.set_xlim(-1,1)
                continue
            if pd.isna( cell ) | (cell in ["plate1-D10","plate1-A7"]) | (cell not in pt_reads.index):
                [_ax.spines[j].set_visible( False ) for j in _ax.spines]
                _ax.set_yticks([])
                _ax.set_xticks([])
                _ax.set_ylim(-1,1)
                _ax.set_xlim(-1,1)
                continue
            values = pt_reads.loc[cell].values
            if values.sum() == 0:
                [_ax.spines[j].set_visible( False ) for j in _ax.spines]
                _ax.set_yticks([])
                _ax.set_xticks([])
                _ax.set_ylim(-1,1)
                _ax.set_xlim(-1,1)
                continue
            _ax.pie( values, colors=cmap )
            _ax.set_ylim(-1,1)
            _ax.set_xlim(-1,1)
        
    for _ax, _sample in zip( _axes[0], all_samples ):
        _ax.set_title( _sample[2:4], fontsize=10, va="center" )
    for row, st in zip( _axes, sample_type_order ):
        row[0].set_ylabel( sample_type_display_dict[st], fontsize=12, rotation="horizontal", ha="right", va="center"  )
    
    legend = [
        Line2D([0], [0], linestyle='none', marker='s', color=color, markeredgecolor="black", markeredgewidth=1, label=label, markersize=10 ) for color, label in zip( cmap, pt_reads.columns )
    ]
    legend = _axes[3,-1].legend( handles=legend, handletextpad=0, loc="center right", frameon=False, fontsize=8, bbox_to_anchor = (0,0,6.5,1), prop={'style': 'italic', 'size':8} )

    _fig.align_titles()

    plt.tight_layout(h_pad=0.0, w_pad=0.1)
    plt.savefig("analysis/plots/figureSX-metagenomic-profile.pdf", bbox_inches="tight" )
    plt.show()
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
