import marimo

__generated_with = "0.23.16"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo
    import pandas as pd
    import matplotlib.pyplot as plt
    from utils import setup_plotting_standards, basic_formatting
    import numpy as np
    from matplotlib.lines import Line2D

    setup_plotting_standards()

    sample_type_order = ['QI', 'BI', 'FP_IS_TCBS', 'RDT_IS_TCBS', 'FP_WS_APW', 'RDT_WS_APW', 'WS_Q', 'WS_SW', 'WS_FP', 'WS_RDT']
    enriched_types = ['QI', 'BI', 'FP_IS_TCBS', 'RDT_IS_TCBS', 'FP_WS_APW', 'RDT_WS_APW']
    sample_type_long_dict = {
        'BI': 'Boiled Isolate',
        'FP_IS_TCBS': 'FP Isolates from TCBS',
        'FP_WS_APW': 'FP WS in APW',
        'QI': 'Qiagen Isolates',
        'RDT_IS_TCBS': 'RDT of Isolates from TCBS',
        'RDT_WS_APW': 'RDT WS in APW',
        'WS_FP': 'WS FP',
        'WS_Q': 'WS Qiagen',
        'WS_RDT': 'WS RDT',
        'WS_SW': 'WS Swabs'
    }
    sample_type_display_dict = {
        'BI': 'Boiled\nIsolate',
        'FP_IS_TCBS': 'FP\nIsolates',
        'FP_WS_APW': 'FP\nAPW',
        'QI': 'Qiagen\nIsolates',
        'RDT_IS_TCBS': 'RDT\nIsolates',
        'RDT_WS_APW': 'RDT\nAPW',
        'WS_FP': 'FP\nWS',
        'WS_Q': 'Qiagen\nWS',
        'WS_RDT': 'RDT\nWS',
        'WS_SW': 'Swab\nWS'
    }
    return (
        Line2D,
        basic_formatting,
        np,
        pd,
        plt,
        sample_type_long_dict,
        sample_type_order,
    )


@app.cell
def _(pd):
    samples = pd.read_csv( "resources/samples.csv" )
    samples.head()
    return (samples,)


@app.cell
def _(pd):
    depth = pd.read_csv( "results/depth.csv.gz",index_col=["chrom","bin"] )
    depth.head()
    return (depth,)


@app.cell
def _(depth, pd, sample_type_order, samples):
    st_depths = list()
    for _st in sample_type_order:
        st_depth = depth[[i for i in samples.loc[samples["sample_type"]==_st,"cell"].to_list() if i in depth.columns]].mean(axis=1)
        st_depth.name = _st
        st_depths.append( st_depth )
    st_depths = pd.concat( st_depths, ignore_index=False, axis=1 )

    #mins = st_depths.min()
    #min_add = np.pow(10, np.floor( np.log10( mins.loc[mins>0].min() ) ) )
    st_depths = st_depths.clip( lower=1e-9 )

    st_depths.head()
    return (st_depths,)


@app.cell
def _(
    Line2D,
    basic_formatting,
    np,
    plt,
    sample_type_long_dict,
    sample_type_order,
    st_depths,
):
    fig, ax = plt.subplots( figsize=(12,2.5), ncols=2, gridspec_kw={"width_ratios" : [2,1]}, sharey=True )
    legend = list()
    for idx, _st in enumerate( sample_type_order ):
        ax[0].plot( _st, data=st_depths.loc["AE003852"], color=f"C{idx}", linewidth=1, zorder=50 )
        ax[1].plot( _st, data=st_depths.loc["AE003853"], color=f"C{idx}", linewidth=1, zorder=50 )
        legend.append( Line2D([0], [0], linestyle='none', marker='s', color=f"C{idx}", markeredgecolor="black", markeredgewidth=1, label=sample_type_long_dict[_st], markersize=8 ) )

    ax[0].set_yscale( "log" )
    ax[0].set_xticks( np.arange(0, 350, 50), labels=[0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0] )
    ax[1].set_xticks( np.arange(0, 125, 25), labels=[0, 0.25, 0.5, 0.75, 1.0] )
    ax[0].set_yticks( [1e-5, 1e-6, 1e-7, 1e-8, 1e-9], labels=["-5", "-6", "-7", "-8", "None"] )
    ax[0].axhline( 1 / 4033501, linewidth=0.5, color='black', linestyle="dashed", zorder=10 )
    ax[1].axhline( 1 / 4033501, linewidth=0.5, color='black', linestyle="dashed", zorder=10 )

    basic_formatting( ax=ax[0], which="both", xlabel="Genome position (Mbp)", ylabel="log Relative coverage", xsize=12, ysize=12, xlims=(-5,305) )
    basic_formatting( ax=ax[1], which="both", xlabel="Genome position (Mbp)", xsize=12, ysize=12, xlims=(-5,112) )

    ax[0].set_title( "Chromosome 1", fontsize=12, fontweight="bold", loc="left" )
    ax[1].set_title( "Chromosome 2", fontsize=12, fontweight="bold", loc="left" )
    ax[1].label_outer()
    ax[1].legend( title="Sample type", handles=legend, handletextpad=0, loc="center right", frameon=False, fontsize=8, alignment="left", title_fontproperties={ "size" : 8, "weight" : "bold" }, bbox_to_anchor = (0,0,1.55,1) )

    plt.tight_layout( w_pad=0.5)
    plt.savefig( "analysis/plots/figureX-depth-comparison-all.pdf" )
    plt.show()
    return


@app.cell
def _(depth):
    depth.loc["AE003852"].tail()
    return


@app.cell
def _(depth):
    one_depth = depth.reset_index()
    one_depth.loc[one_depth["chrom"]=="AE003853","bin"] += one_depth.loc[one_depth["chrom"]=="AE003852","bin"].max() + 1
    one_depth = (one_depth
        .set_index( "bin" )
        .drop( columns=["chrom"] )
        .clip( lower=1e-9, upper=1e-5 )
        .dropna( axis=1 ))
    return (one_depth,)


@app.cell
def _(
    basic_formatting,
    np,
    one_depth,
    plt,
    sample_type_long_dict,
    sample_type_order,
    samples,
):
    _fig, _axes = plt.subplots( figsize=(12,10), ncols=2, nrows=5, sharex=True, sharey=True )
    for _st, _ax in zip( sample_type_order, _axes.flatten() ):
        st_cells = samples.loc[samples["sample_type"]==_st,"cell"].to_list()
        for _cell in st_cells:
            if _cell not in one_depth.columns:
                continue
            _ax.plot( _cell, data=one_depth, color='black', linewidth=1, alpha=0.1, zorder=100 )
        _ax.axhline( 1 / 4033501, linewidth=0.75, color='red', linestyle="dashed", zorder=150 )
        _ax.axvline( 297, linewidth=1, color="black", zorder=200 )
        _ax.set_yscale( "log" )
        _ax.set_xticks( np.arange(0, 450, 50), labels=[f"{i/100:.1f}" for i in np.arange(0, 450, 50)] )
        _ax.set_yticks( [1e-5, 1e-6, 1e-7, 1e-8, 1e-9], labels=["-5", "-6", "-7", "-8", "None"] )
        basic_formatting( _ax, xlabel="Genomic position (Mbp)", ylabel="log Relative coverage", which="both", xlims=(-15, 415))
        _ax.set_title( sample_type_long_dict[_st], fontsize=12, fontweight="bold", loc="left" )
        _ax.label_outer()
    plt.tight_layout()
    plt.savefig( "analysis/plots/figureX-depth-comparison-individual.pdf")
    plt.show()
    return


@app.cell
def _(pd):
    dep = pd.read_csv( "resources/plate3-B9.depth", sep="\t", header=None, names=["chromosome", "pos", "coverage"], index_col=[0,1] )
    return (dep,)


@app.cell
def _(dep):
    dep.loc["AE003852"].plot( xlim=(401_000,407_500), figsize=(10,4) )
    return


@app.cell
def _(dep):
    dep.loc["AE003853"].plot( xlim=(528_000,545_000), figsize=(10,4) )
    return


if __name__ == "__main__":
    app.run()
