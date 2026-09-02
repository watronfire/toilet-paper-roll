import marimo

__generated_with = "0.23.16"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo
    import marimo as mo
    import pandas as pd
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    import matplotlib.ticker as mticker
    from utils import setup_plotting_standards, basic_formatting, get_sample_types
    import numpy as np
    from statsmodels.stats.proportion import proportion_confint
    import seaborn as sns

    setup_plotting_standards()

    sample_type_order, sample_type_long_dict, sample_type_display_dict = get_sample_types()

    samples = pd.read_csv( "resources/samples.csv" )
    return (
        basic_formatting,
        mticker,
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
    var = pd.read_csv( "results/variants_subsample.csv", header=None, names=["cell", "reads", "pi", "pi_median", "af", "af_median", "variant_count"] )
    var = var.merge( samples, on="cell" )
    var["reads"] = var["reads"] / 2
    # This sample is probably a non-cholerae vibrio. Should be improved by filtering approach but isn't.
    var = var.loc[var["cell"]!="plate1-H7"]
    var
    return (var,)


@app.cell
def _(
    basic_formatting,
    mticker,
    np,
    plt,
    sample_type_display_dict,
    sample_type_order,
    sns,
    var,
):
    # 1. Common boxplot styling
    _cmap = {"Isolates" : "#66c2a5", "APW" : "#fc8d62", "Stool" : "#8da0cb"}
    _style_kwargs = {
        #"color": "#E5E4E2",
        "linecolor": "black",
        "fliersize": 5,
        "linewidth": 0.75,
        "formatter": lambda x: sample_type_display_dict[x].split("\n")[0],
        "medianprops": {"linewidth": 1},
        "capwidths": 0.15,
    }

    # 2. Define target variables and their panel-specific settings per row
    _row_configs = [
        {
            "y_var": "variant_count",
            "ylabel": "Within-host variant counts",
            "ylims": (-1, 40.5),
            "yticks_major": range(0, 50, 10),
            "yticks_minor": range(0, 40, 2),
        },
        {
            "y_var": "af",
            "ylabel": "Within-host variant frequency",
            "ylims": (-0.01, 0.51),  # Adjust based on your frequency range (e.g., 0 to 1)
            "yticks_major": np.arange(0,0.6,0.1),
            "yticks_minor": np.arange(0,0.5,0.02),
        },
    ]

    # Column subsets and x-labels
    _col_config = [
        (0, sample_type_order[:4], "Isolates"),
        (1, sample_type_order[4:6], "APW"),
        (2, sample_type_order[6:], "Stool"),
    ]

    # 3. Create a 2x3 subplot grid (sharing y-axis per row)
    _fig, _axes = plt.subplots(
        nrows=2,
        ncols=3,
        figsize=(8, 7),
        gridspec_kw={"width_ratios": [4, 2, 4]}
    )

    # 4. Iterate over rows and columns
    for _row_idx, _r_cfg in enumerate(_row_configs):
        for _col_idx, _sample_subset, _x_label in _col_config:
            _ax = _axes[_row_idx, _col_idx]
            _subset_data = var[var["sample_type"].isin(_sample_subset)]

            # Plot boxplot
            sns.boxplot(
                data=_subset_data,
                x="sample_type",
                y=_r_cfg["y_var"],
                order=_sample_subset,
                ax=_ax,
                color=_cmap[_x_label],
                flierprops={"color": "#E5E4E2", "markeredgecolor": "black", "markeredgewidth": 0.75, "markerfacecolor": _cmap[_x_label]},
                **_style_kwargs,
            )

            # Apply formatting (X label on bottom row only, Y label on left col via label_outer)
            basic_formatting(
                _ax,
                ylabel=_r_cfg["ylabel"],
                xlabel=_x_label,
                ylims=_r_cfg["ylims"],
                xsize=12,
                ysize=12,
                spines=[],
            )
            if _row_idx == 1:
                _ax.yaxis.set_major_formatter( mticker.PercentFormatter(1) )
            _ax.tick_params(axis="both", which="both", size=0)
            if _col_idx > 0: 
                _ax.set_yticklabels([])
                _ax.set_ylabel( "" )

            #if _col_idx == 1:
            #    _ax.axhspan(-1,40, color="black", alpha=0.04, zorder=2, linewidth=0)

            _ax.set_yticks(_r_cfg["yticks_major"])
            _ax.set_yticks(_r_cfg["yticks_minor"], minor=True)
            _ax.set_xlabel( _x_label, fontweight="medium" )


    plt.tight_layout(w_pad=0.5)
    plt.savefig( "analysis/plots/figure3-within-host-variants.pdf" )
    plt.show()
    return


if __name__ == "__main__":
    app.run()
