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

    setup_plotting_standards()

    sample_type_order, sample_type_long_dict, sample_type_display_dict = get_sample_types()

    samples = pd.read_csv( "resources/samples.csv" )
    return (
        Line2D,
        basic_formatting,
        mticker,
        np,
        pd,
        plt,
        sample_type_display_dict,
        sample_type_order,
        samples,
    )


@app.cell
def _(pd, samples):
    krak = pd.read_csv( "results/reports/kraken2_reports.csv" )
    krak = krak.merge( samples, on="cell", how="left" )

    krak["total_reads"] = krak["classified"] + krak["unclassified"]
    krak["fraction_unclassified"] = krak["unclassified"] / krak["total_reads"]
    krak["fraction_human"] = krak["S_homo-sapiens"] / krak["total_reads"]
    krak["fraction_vibrio"] = krak["O_vibrionales"] / krak["total_reads"]
    #krak = krak.loc[~((krak["sample_type"]=="WS_RDT")&(krak["sequencing_run"]=="Filter_Paper_Sequencing_Run_1_061325"))]
    krak = krak.groupby( ["sample", "sample_type"] ).first().reset_index()

    krak_samples = krak["sample"].sort_values().unique()
    krak
    return krak, krak_samples


@app.cell
def _(
    Line2D,
    basic_formatting,
    krak,
    krak_samples,
    mticker,
    np,
    plt,
    sample_type_display_dict,
    sample_type_order,
):
    _fig, _axes = plt.subplots( dpi=200, figsize=(5,7), nrows=3, height_ratios=[4,2,4], sharex=True )

    cmap={
        "fraction_vibrio" : "#66c4ac",
        "fraction_human" : "#9ed1f1",
        "fraction_unclassified" : "gainsboro", 
        "Other" : "#E69F00" 
    } 

    _col_config = [
        (0, sample_type_order[:4], "Isolates"),
        (1, sample_type_order[4:6], "APW"),
        (2, sample_type_order[6:], "Stool"),
    ]

    _width = 1/(len(krak_samples)+2)
    for _ax, (_col_idx, _sample_subset, _x_label) in zip( _axes, _col_config ):
        for _idx, _typ in enumerate( _sample_subset ):
            for _jdx, _sample in enumerate( krak_samples ):
                _position = _idx + (_jdx / (len(krak_samples) + 2))
                _data = krak.loc[(krak["sample"]==_sample)&(krak["sample_type"]==_typ)]
                _bottom = 0
                for _column in ["fraction_vibrio", "fraction_human", "fraction_unclassified"]:
                    _ax.barh( _position, _data[_column], left=_bottom, height=_width, color=cmap[_column], zorder=10, align="edge" )
                    _bottom += _data[_column]
                _ax.barh( _position, 1-_bottom, left=_bottom, color=cmap["Other"], height=_width, zorder=10, align="edge")
        
        for _idx in range( 1, len( _sample_subset ) ):
            _ax.axhline( _idx-0.04, color="black", linewidth=1, zorder=100 )
    
        _ax.set_xticks(np.arange(0, 1.2, 0.2) )
        _ax.set_xticks(np.arange(0, 1.05, 0.05), minor=True )
        _ax.set_yticks( np.arange( 0.5, len( _sample_subset )+0.5, 1 ), [sample_type_display_dict[i].split( "\n" )[0] for i in _sample_subset] )
    
        basic_formatting( _ax, which="x", xlims=(-0.01,1.02), ylims=(-0.1,len(_sample_subset)), xsize=12, ysize=12, spines=[] )
        _ax.set_ylabel( _x_label, fontweight="medium", )
        _ax.label_outer()
        _ax.xaxis.set_major_formatter( mticker.PercentFormatter(1) )
        _ax.yaxis.set_inverted(True)

    _legend_props = {"linestyle" : "none", "marker" : "s", "markeredgecolor" : "black", "markeredgewidth" : 1, "markersize" : 10}

    _legend = [
        Line2D([0], [0], color=cmap["fraction_vibrio"], label="Vibrio", **_legend_props ),
        Line2D([0], [0], color=cmap["fraction_human"], label="Human", **_legend_props ),
        Line2D([0], [0], color=cmap["fraction_unclassified"], label="Unclassified", **_legend_props ),
        Line2D([0], [0], color=cmap["Other"], label="Other", **_legend_props ),
    ]

    _legend = _axes[0].legend( handles=_legend, handletextpad=0, loc="upper center", frameon=False, fontsize=10, bbox_to_anchor = (0,0,1,1.15), ncols=4 )

    _fig.align_ylabels()
    plt.tight_layout(h_pad=0.3)
    plt.savefig( "analysis/plots/figure4-metagenomic-profile.pdf" )
    plt.show()
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
