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
    from scipy.special import logit, expit

    import pytensor
    pytensor.config.cxx = ""
    import pymc as pm
    import arviz as az

    setup_plotting_standards()
    plt.rcParams['mathtext.fontset'] = 'custom'

    sample_type_order, sample_type_long_dict, sample_type_display_dict = get_sample_types()

    sample_type_display_dict = {
        'BI': 'Isolate - Boiled',
        'FP_IS_TCBS': 'Isolate - Filter paper',
        'FP_WS_APW': 'APW - Filter paper',
        'QI': 'Isolate - Qiagen',
        'RDT_IS_TCBS': 'Isolate - RDT',
        'RDT_WS_APW': 'APW - RDT',
        'WS_FP': 'Stool - Filter paper',
        'WS_Q': 'Stool - Qiagen',
        'WS_RDT': 'Stool - RDT',
        'WS_SW': 'Stool - Swab'
    }

    st_numeric = {st : i for i, st in enumerate(sample_type_order)}

    samples = pd.read_csv( "resources/samples.csv" )
    return (
        az,
        basic_formatting,
        mticker,
        np,
        pd,
        plt,
        pm,
        sample_type_display_dict,
        sample_type_order,
        samples,
        st_numeric,
    )


@app.cell
def _(pd):
    ct = pd.read_csv( "resources/WGS-qPCR.csv", usecols=["cell", "O1", "toxR", "ctxA"] )
    ct = ct.dropna( subset=["cell"] )
    ct.head()
    return (ct,)


@app.cell
def _(ct, np, pd, samples, st_numeric):
    rare = pd.read_csv( "results/coverage_subsample.csv", header=None, names=["cell","reads", "coverage", "median_depth"])
    rare = rare.merge( samples, on="cell", how="left" ).merge( ct, on="cell", how="left" )
    rare = rare.loc[rare["reads"].isin([100_000,250_000,500_000,1_000_000,2_500_000,5_000_000])].dropna( subset=["O1", "coverage"] )
    rare["sample_type_numeric"] = rare["sample_type"].apply( lambda x: st_numeric[x] )
    rare["coverage"] = rare["coverage"].clip(1e-6, 1 - 1e-6)
    rare["log_reads"] = np.log( rare["reads"] )
    rare
    return (rare,)


@app.cell
def _(pm, rare, sample_type_order):
    with pm.Model(coords={"sample_type": sample_type_order}) as rarefaction_model:

        # --- Data ---
        _ct_values       = pm.Data("ct", rare["O1"].values)
        _log_reads       = pm.Data("log_reads", rare["log_reads"] )  # work on log scale
        _coverage        = pm.Data("coverage", rare["coverage"].values)
        _sample_type_idx = pm.Data("sample_type_idx", rare["sample_type_numeric"].values)

        # --- Hyperpriors for alpha (log-reads midpoint intercept) ---
        # log(1e6) ≈ 13.8: prior centered on ~1M reads for 50% coverage at Ct=0
        # the Ct term (beta * ct) will shift this to realistic values
        _mu_alpha    = pm.Normal("mu_alpha", mu=13, sigma=3)
        _sigma_alpha = pm.HalfNormal("sigma_alpha", sigma=2)

        # --- Hyperpriors for beta (how midpoint shifts per Ct unit) ---
        # Positive: higher Ct → midpoint shifts right → more reads needed
        _mu_beta    = pm.Normal("mu_beta", mu=0.3, sigma=0.5)
        _sigma_beta = pm.HalfNormal("sigma_beta", sigma=0.5)

        # --- Hyperpriors for k_r (rarefaction curve steepness) ---
        # Positive, so use HalfNormal; steepness on log-read scale
        # kr on log scale to prevent sign flips
        _mu_log_kr    = pm.Normal("mu_log_kr",    mu=0,  sigma=1)
        _sigma_log_kr = pm.HalfNormal("sigma_log_kr",    sigma=0.5)

        # --- Partially pooled per-sample-type parameters ---
        _alpha_offset = pm.Normal("alpha_offset", mu=0, sigma=1, dims="sample_type")
        _alpha = pm.Deterministic("alpha", _mu_alpha + _alpha_offset * _sigma_alpha, dims="sample_type")

        _beta_offset = pm.Normal("beta_offset", mu=0, sigma=1, dims="sample_type")
        _beta = pm.Deterministic("beta", _mu_beta + _beta_offset * _sigma_beta, dims="sample_type")

        _kr_offset = pm.Normal("kr_offset",    mu=0, sigma=1, dims="sample_type")
        _kr = pm.Deterministic("kr", pm.math.exp(_mu_log_kr + _kr_offset * _sigma_log_kr), dims="sample_type")

        # --- Midpoint in log-read space as a function of Ct ---
        _midpoint = _alpha[_sample_type_idx] + _beta[_sample_type_idx] * _ct_values

        # --- Logistic in log(reads) ---
        _mu = pm.Deterministic( "mu", 1 / (1 + pm.math.exp(-_kr[_sample_type_idx] * (_log_reads - _midpoint) ) ) )

        # --- Beta likelihood: keeps predictions in (0, 1) ---
        _kappa = pm.HalfNormal("kappa", sigma=1)
        _y = pm.Normal( "y", mu=_mu, sigma=_kappa, observed=_coverage)

        rarefaction_trace = pm.sample(draws=1000, tune=5000, return_inferencedata=True, nuts_sampler="nutpie" )
        rarefaction_trace = pm.compute_log_likelihood(rarefaction_trace)
        rarefaction_trace = pm.sample_posterior_predictive(rarefaction_trace, sample_vars=["y", "mu"], extend_inferencedata=True)
    return (rarefaction_trace,)


@app.cell
def _(az, rarefaction_trace):
    az.summary( rarefaction_trace ).sort_values( "r_hat", ascending=False )
    return


@app.cell
def _(
    mticker,
    np,
    plt,
    rarefaction_trace,
    sample_type_display_dict,
    sample_type_order,
):
    ct_grid    = np.linspace(0, 40, 50)
    reads_grid = np.logspace(4, 9, 50)
    log_reads_grid = np.log(reads_grid)

    X, Y = np.meshgrid( ct_grid, log_reads_grid )
    data_extent = [ct_grid.min(), ct_grid.max(), log_reads_grid.min(), log_reads_grid.max()]

    log_ticks = np.log( np.array([d * 10**e for e in range(4, 9) for d in range(1, 10)]) )

    _fig, _axes = plt.subplots( dpi=200, figsize=(11,4), nrows=2, ncols=5, sharex=True, sharey=True )

    for _ax, _st_name in zip( _axes.flatten(), sample_type_order ):
        _alpha_post = rarefaction_trace.posterior["alpha"].sel(sample_type=_st_name).stack(sample=("chain", "draw")).values   # (chain, draw)
        _beta_post  = rarefaction_trace.posterior["beta"].sel(sample_type=_st_name).stack(sample=("chain", "draw")).values 
        _kr_post    = rarefaction_trace.posterior["kr"].sel(sample_type=_st_name).stack(sample=("chain", "draw")).values 

        # midpoint[draw, ct]: shape (draws, n_ct)
        _midpoint = _alpha_post[:, None] + _beta_post[:, None] * ct_grid[None, :]

        # pred: (n_reads, n_samples, n_ct)
        _pred = 1 / (1 + np.exp(-_kr_post[None, :, None] * (log_reads_grid[:, None, None] - _midpoint[None, :, :])))

        # posterior mean coverage surface [reads x ct]
        _pred_mean = _pred.mean(axis=1)

        _im = _ax.imshow( _pred_mean, aspect="auto", cmap="YlGnBu", origin="lower", extent=data_extent, vmin=0, vmax=1 )
        _ax.set_ylim(np.log( 5e4 ),np.log( 4e8 ) )
        _ax.set_xlim(0,40 )

        _CS = _ax.contour( X, Y, _pred_mean, levels=[0.5, 0.95], colors=["white", "white"], extent=data_extent, linewidths=1, linestyles=["dashed","solid"] )
        _ax.clabel(_CS, _CS.levels, fmt=lambda x: f"{x:.0%}", fontsize=10, inline_spacing=6, inline=True )

        _ax.set_title( sample_type_display_dict[_st_name], fontsize=9, fontweight="medium", loc="left", pad=5 )
        _ax.set_yticks( np.log([1e5, 1e6,1e7, 1e8]), ["$10^5$", "$10^6$", "$10^7$", "$10^8$" ] )
        _ax.set_yticks( log_ticks, minor=True )
        _ax.set_xticks( np.arange(0,45,10 ) )
        _ax.set_xticks( np.arange(0,40,1 ), minor=True )
        _ax.set_ylabel( "Reads", fontsize=10 )
        _ax.set_xlabel( "O1 Ct value", fontsize=10 )
        _ax.label_outer()
        _ax.tick_params( axis="both", which="both", labelsize=10, color="silver" )
        _ax.tick_params( axis="both", which="major", length=5 )
        _ax.set_ylim(np.log( 3e4 ),np.log( 4e8 ) )
        [_ax.spines[j].set_visible( False ) for j in _ax.spines]

    _cbar = _fig.colorbar(_im, ax=_axes.ravel(), shrink=0.5, orientation='vertical', format=mticker.PercentFormatter(1), fraction=0.05, pad=0.03, aspect=12 )
    _cbar.set_label('Coverage', rotation=270, fontsize=10, labelpad=10 )
    _cbar.outline.set_visible(False)
    _cbar.ax.tick_params( axis="both", direction="in", left=True, right=True, labelsize=10, color="#E5E4E2" )

    #plt.tight_layout()
    _fig.savefig( "analysis/plots/figure5-coverage-vs-reads-vs-ct.pdf", bbox_inches="tight" )
    plt.show()
    return X, Y, ct_grid, data_extent, log_reads_grid, log_ticks


@app.cell
def _(
    X,
    Y,
    ct_grid,
    data_extent,
    log_reads_grid,
    log_ticks,
    mticker,
    np,
    plt,
    rarefaction_trace,
    sample_type_display_dict,
    sample_type_order,
):
    _fig, _axes = plt.subplots( dpi=200, figsize=(11,4), nrows=2, ncols=5, sharex=True, sharey=True )

    for _ax, _st_name in zip( _axes.flatten(), sample_type_order ):
        _alpha_post = rarefaction_trace.posterior["alpha"].sel(sample_type=_st_name).stack(sample=("chain", "draw")).values   # (chain, draw)
        _beta_post  = rarefaction_trace.posterior["beta"].sel(sample_type=_st_name).stack(sample=("chain", "draw")).values 
        _kr_post    = rarefaction_trace.posterior["kr"].sel(sample_type=_st_name).stack(sample=("chain", "draw")).values 

        # midpoint[draw, ct]: shape (draws, n_ct)
        _midpoint = _alpha_post[:, None] + _beta_post[:, None] * ct_grid[None, :]

        # pred: (n_reads, n_samples, n_ct)
        _pred = 1 / (1 + np.exp(-_kr_post[None, :, None] * (log_reads_grid[:, None, None] - _midpoint[None, :, :])))

        # posterior mean coverage surface [reads x ct]
        _succ = (_pred > 0.9).mean(axis=1)

        _im = _ax.imshow( _succ, aspect="auto", cmap="YlGnBu", origin="lower", extent=data_extent, vmin=0, vmax=1 )
        _ax.set_ylim(np.log( 5e4 ),np.log( 4e8 ) )
        _ax.set_xlim(0,40 )

        _CS = _ax.contour( X, Y, _succ, levels=[0.5, 0.95], colors=["white", "white"], extent=data_extent, linewidths=1, linestyles=["dashed","solid"] )
        _ax.clabel(_CS, _CS.levels, fmt=lambda x: f"{x:.0%}", fontsize=10, inline_spacing=6, inline=True )

        _ax.set_title( sample_type_display_dict[_st_name], fontsize=9, fontweight="medium", loc="left", pad=5 )
        _ax.set_yticks( np.log([1e5, 1e6,1e7, 1e8]), ["$10^5$", "$10^6$", "$10^7$", "$10^8$" ] )
        _ax.set_yticks( log_ticks, minor=True )
        _ax.set_xticks( np.arange(0,45,10 ) )
        _ax.set_xticks( np.arange(0,40,1 ), minor=True )
        _ax.set_ylabel( "Reads", fontsize=10 )
        _ax.set_xlabel( "O1 Ct value", fontsize=10 )
        _ax.label_outer()
        _ax.tick_params( axis="both", which="both", labelsize=10, color="silver" )
        _ax.tick_params( axis="both", which="major", length=5 )
        _ax.set_ylim(np.log( 3e4 ),np.log( 4e8 ) )
        [_ax.spines[j].set_visible( False ) for j in _ax.spines]
        #_ax.grid( axis="both", which="major", color="#CDCDCD" )

    _cbar = _fig.colorbar(_im, ax=_axes.ravel(), shrink=0.5, orientation='vertical', format=mticker.PercentFormatter(1), fraction=0.05, pad=0.03, aspect=12 )
    _cbar.set_label('Pr(complete genome)', rotation=270, fontsize=10, labelpad=10 )
    _cbar.outline.set_visible(False)
    _cbar.ax.tick_params( axis="both", direction="in", left=True, right=True, labelsize=10, color="#E5E4E2" )

    #plt.tight_layout()
    _fig.savefig( "analysis/plots/figure5-sequencing-success-vs-reads-vs-ct.pdf", bbox_inches="tight" )
    plt.show()
    return


@app.cell
def _(
    az,
    ct_grid,
    mticker,
    np,
    plt,
    rare,
    rarefaction_trace,
    sample_type_display_dict,
    sample_type_order,
):
    for _mod in [1, 5]:
        _fig, _axes = plt.subplots( dpi=200, figsize=(11,4.5), nrows=2, ncols=5, sharex=True, sharey=True )
        _reads = 1_000_000 * _mod
    
        for _ax, _st_name in zip( _axes.flatten(), sample_type_order ):
            _alpha_post = rarefaction_trace.posterior["alpha"].sel(sample_type=_st_name).stack(sample=("chain", "draw")).values   # (chain, draw)
            _beta_post  = rarefaction_trace.posterior["beta"].sel(sample_type=_st_name).stack(sample=("chain", "draw")).values 
            _kr_post    = rarefaction_trace.posterior["kr"].sel(sample_type=_st_name).stack(sample=("chain", "draw")).values 
    
            # midpoint[draw, ct]: shape (draws, n_ct)
            _midpoint = _alpha_post[:, None] + _beta_post[:, None] * ct_grid[None, :]
    
            # pred: (n_reads, n_samples, n_ct)
            _pred = 1 / (1 + np.exp(-_kr_post[None, :, None] * (np.log(_reads) - _midpoint[None, :, :])))
            _pred_mean = _pred.mean(axis=1)[0]
            _pred_hdi = az.hdi(_pred[0].T, prob=0.95 )
    
            _ax.plot( ct_grid, _pred_mean, color='black', linewidth=1, linestyle='dashed', zorder=10 )
            _ax.fill_between( ct_grid, _pred_hdi[:,0], _pred_hdi[:,1], color='black', linewidth=0, alpha=0.1, zorder=5 )
    
            _points = _ax.scatter( "O1", "coverage", data=rare.loc[(rare["sample_type"]==_st_name)&(rare["reads"]==_reads)], zorder=20, alpha=0.4, linewidth=0)
            _points.set_clip_on( False )
            _ax.set_title( sample_type_display_dict[_st_name], fontsize=8, fontweight="bold", loc="left" )
            _ax.tick_params( axis="both", which="both", labelsize=10, length=0)
            _ax.set_yticks( np.arange(0,1.25,0.25))
            _ax.set_yticks( np.arange(0,1.0,0.05), minor=True)
            _ax.set_xticks( np.arange(0,45,10))
            _ax.set_xticks( np.arange(0,40,2), minor=True)
            _ax.set_ylim(-0.01,1.01)
            _ax.set_xlim(-0.1,40.1)
    
            [_ax.spines[j].set_visible( False ) for j in _ax.spines]
            _ax.yaxis.set_major_formatter( mticker.PercentFormatter(1) )
            _ax.set_xlabel( "O1 Ct value", fontsize=10)
            _ax.set_ylabel( f"Coverage w/\n{_mod} million reads", fontsize=10)
            _ax.grid( which="major", color="#efefef", linewidth=1 )
            _ax.grid( which="minor", color="#efefef", linewidth=0.5 )
            _ax.label_outer()
    
        plt.tight_layout()
        plt.savefig( f"analysis/plots/figureSX-coverage-vs-ct-{_mod}mil.pdf")
        plt.show()
    return


@app.cell
def _(az, rarefaction_trace):
    az.residual_r2( data=rarefaction_trace, pred_mean="mu" )
    return


@app.cell
def _(basic_formatting, np, plt, rare, rarefaction_trace):
    # Posterior mean prediction for each observation
    mu_samples = (1 / (1 + np.exp(
        -rarefaction_trace.posterior["kr"].values[:, :, rare["sample_type_numeric"].values] *
        (rare["log_reads"].values -
         rarefaction_trace.posterior["alpha"].values[:, :, rare["sample_type_numeric"].values] -
         rarefaction_trace.posterior["beta"].values[:, :, rare["sample_type_numeric"].values] *
         rare["O1"].values)
    ))).mean(axis=(0, 1))  # average over chains and draws

    residuals = rare["coverage"].values.clip(1e-6, 1 - 1e-6) - mu_samples

    fig, axes = plt.subplots(1, 3, figsize=(12, 4))
    axes[0].scatter(rare["O1"],  residuals, alpha=0.4, linewidth=0, zorder=100)
    axes[1].scatter(rare["log_reads"], residuals, alpha=0.4, linewidth=0, zorder=100)
    axes[2].scatter(mu_samples, residuals, alpha=0.4, linewidth=0, zorder=100)
    for ax, _label in zip( axes, ["O1 Ct value", "log Reads", "Fitted value"] ):
        ax.axhline(0, color="red", linestyle="--", zorder=150)
        ax.label_outer()
        basic_formatting( ax, spines=[], which="both", xlabel=_label, ylabel="Residual")

    plt.tight_layout()
    plt.savefig( "analysis/plots/figureSX-coverage-vs-ct-residuals.pdf" )
    plt.show()
    return


if __name__ == "__main__":
    app.run()
