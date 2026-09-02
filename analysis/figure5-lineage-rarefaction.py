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
    from sklearn.calibration import calibration_curve

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
        calibration_curve,
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
    vibe = pd.read_csv( "results/lineage_subsample.csv" )
    vibe = vibe.merge( samples, on="cell", how="left" ).merge( ct, on="cell", how="left" )
    vibe["sample_type_numeric"] = vibe["sample_type"].map( st_numeric )
    vibe = vibe.loc[vibe["reads"].isin([100_000,250_000,500_000,1_000_000,2_500_000,5_000_000])].dropna( subset=["O1", "correct"] )
    vibe["log_reads"] = np.log( vibe["reads"] )
    vibe = vibe.loc[vibe["O1"]<40] # Turns out I can't conservatively estimate these as >40 Ct values.
    vibe.head()
    return (vibe,)


@app.cell
def _(pm, sample_type_order, vibe):
    with pm.Model(coords={"sample_type": sample_type_order}) as vibes_model:

        # --- Data ---
        _ct_values       = pm.Data("ct", vibe["O1"].values)
        _log_reads       = pm.Data("log_reads", vibe["log_reads"].values)
        _correct         = pm.Data("correct", vibe["correct"].values)
        _sample_type_idx = pm.Data("sample_type_idx", vibe["sample_type_numeric"].values)

        # --- Hyperpriors ---
        _mu_alpha    = pm.Normal("mu_alpha", mu=13, sigma=3)
        _sigma_alpha = pm.HalfNormal("sigma_alpha", sigma=2)

        # mu_beta is now strictly a POSITIVE multiplier representing 
        # "how many additional log-reads are needed per 1 unit increase in Ct"
        _mu_beta     = pm.HalfNormal("mu_beta", sigma=0.5) 
        _sigma_beta  = pm.HalfNormal("sigma_beta", sigma=0.5)

        _mu_log_kr    = pm.Normal("mu_log_kr", mu=0, sigma=1)
        _sigma_log_kr = pm.HalfNormal("sigma_log_kr", sigma=0.5)

        # --- Partially Pooled Parameters ---
        _alpha_offset = pm.Normal("alpha_offset", mu=0, sigma=1, dims="sample_type")
        _alpha = pm.Deterministic("alpha", _mu_alpha + _alpha_offset * _sigma_alpha, dims="sample_type")

        # Using Bound/HalfNormal to force beta positive, meaning Ct pushes the midpoint RIGHT
        _beta_offset = pm.Normal("beta_offset", mu=0, sigma=1, dims="sample_type")
        _beta = pm.Deterministic("beta", _mu_beta + _beta_offset * _sigma_beta, dims="sample_type")

        _kr_offset = pm.Normal("kr_offset", mu=0, sigma=1, dims="sample_type")
        _kr = pm.Deterministic("kr", pm.math.exp(_mu_log_kr + _kr_offset * _sigma_log_kr), dims="sample_type")

        # --- The Mathematical Lock-In ---
        # We define the midpoint where higher Ct explicitly INCREASES the required log_reads
        _midpoint = _alpha[_sample_type_idx] + _beta[_sample_type_idx] * _ct_values

        # Standard invlogit: log_reads minus the shifting midpoint, scaled by kr
        _mu_logit = _kr[_sample_type_idx] * (_log_reads - _midpoint)
        _mu = pm.Deterministic( "mu", pm.math.invlogit(_mu_logit) )

        _y = pm.Bernoulli("y", p=_mu, observed=_correct)

        vibes_trace = pm.sample(draws=1000, tune=1000, return_inferencedata=True, nuts_sampler="nutpie" )
        vibes_trace = pm.compute_log_likelihood(vibes_trace)
        vibes_trace = pm.sample_posterior_predictive(vibes_trace, sample_vars=["y", "mu"], extend_inferencedata=True)
    return (vibes_trace,)


@app.cell
def _(az, vibes_trace):
    az.summary( vibes_trace ).sort_values( "r_hat", ascending=False )
    return


@app.cell
def _(
    mticker,
    np,
    plt,
    sample_type_display_dict,
    sample_type_order,
    vibes_trace,
):
    ct_grid    = np.linspace(0, 40, 50)
    reads_grid = np.logspace(3, 8, 50)
    log_reads_grid = np.log(reads_grid)

    X, Y = np.meshgrid( ct_grid, log_reads_grid )
    data_extent = [ct_grid.min(), ct_grid.max(), log_reads_grid.min(), log_reads_grid.max()]

    log_ticks = np.log( np.array([d * 10**e for e in range(3, 9) for d in range(1, 10)]) )

    _fig, _axes = plt.subplots( figsize=(11,4), nrows=2, ncols=5, sharex=True, sharey=True )

    for _ax, _st_name in zip( _axes.flatten(), sample_type_order ):
        _alpha_post = vibes_trace.posterior["alpha"].sel(sample_type=_st_name).stack(sample=("chain", "draw")).values   # (chain, draw)
        _beta_post  = vibes_trace.posterior["beta"].sel(sample_type=_st_name).stack(sample=("chain", "draw")).values 
        _kr_post    = vibes_trace.posterior["kr"].sel(sample_type=_st_name).stack(sample=("chain", "draw")).values 

        # midpoint[draw, ct]: shape (draws, n_ct)
        _midpoint = _alpha_post[:, None] + _beta_post[:, None] * ct_grid[None, :]

        # pred: (n_reads, n_samples, n_ct)
        _pred = 1 / (1 + np.exp(-_kr_post[None, :, None] * (log_reads_grid[:, None, None] - _midpoint[None, :, :])))

        # posterior mean coverage surface [reads x ct]
        _pred_mean = _pred.mean(axis=1)

        _im = _ax.imshow( _pred_mean, aspect="auto", cmap="YlGnBu", origin="lower", extent=data_extent, vmin=0, vmax=1 )
        _ax.set_ylim(np.log( 3e3 ),np.log( 4e7 ) )
        _ax.set_xlim(0,40 )

        _CS = _ax.contour( X, Y, _pred_mean, levels=[0.5, 0.95], colors=["white", "white"], extent=data_extent, linewidths=1, linestyles=["dashed","solid"] )
        _ax.clabel(_CS, _CS.levels, fmt=lambda x: f"{x:.0%}", fontsize=10, inline_spacing=6, inline=True )

        _ax.set_title( sample_type_display_dict[_st_name], fontsize=9, fontweight="medium", loc="left", pad=5 )
        _ax.set_yticks( np.log([1e4,1e5, 1e6,1e7]), ["$10^4$", "$10^5$", "$10^6$", "$10^7$" ] )
        _ax.set_yticks( log_ticks, minor=True )
        _ax.set_xticks( np.arange(0,45,10 ) )
        _ax.set_xticks( np.arange(0,40,1 ), minor=True )
        _ax.set_ylabel( "Reads", fontsize=10 )
        _ax.set_xlabel( "O1 Ct value", fontsize=10 )
        _ax.label_outer()
        _ax.tick_params( axis="both", which="both", labelsize=10, color="silver" )
        _ax.tick_params( axis="both", which="major", length=5 )
        _ax.set_ylim(np.log( 3e3 ),np.log( 4e7 ) )
        #_ax.set_ylim(np.log( 3e4 ),np.log( 4e8 ) )
        [_ax.spines[j].set_visible( False ) for j in _ax.spines]
    
    _cbar = _fig.colorbar(_im, ax=_axes.ravel(), shrink=0.5, orientation='vertical', format=mticker.PercentFormatter(1), fraction=0.05, pad=0.03, aspect=12 )
    _cbar.set_label('Pr(Correct lineage)', rotation=270, fontsize=10, labelpad=10 )
    _cbar.outline.set_visible(False)
    _cbar.ax.tick_params( axis="both", direction="in", left=True, right=True, labelsize=10, color="#E5E4E2" )

    #plt.tight_layout()
    _fig.savefig( "analysis/plots/figure5-lineage-vs-reads-vs-ct.pdf", bbox_inches="tight" )
    plt.show()
    return (ct_grid,)


@app.cell
def _(
    az,
    ct_grid,
    mticker,
    np,
    plt,
    sample_type_display_dict,
    sample_type_order,
    vibe,
    vibes_trace,
):
    for _mod in [1, 5]:
        _fig, _axes = plt.subplots( dpi=200, figsize=(11,4.5), nrows=2, ncols=5, sharex=True, sharey=True )
        _reads = 1_000_000 * _mod
    
        for _ax, _st_name in zip( _axes.flatten(), sample_type_order ):
            _alpha_post = vibes_trace.posterior["alpha"].sel(sample_type=_st_name).stack(sample=("chain", "draw")).values   # (chain, draw)
            _beta_post  = vibes_trace.posterior["beta"].sel(sample_type=_st_name).stack(sample=("chain", "draw")).values 
            _kr_post    = vibes_trace.posterior["kr"].sel(sample_type=_st_name).stack(sample=("chain", "draw")).values 
    
            # midpoint[draw, ct]: shape (draws, n_ct)
            _midpoint = _alpha_post[:, None] + _beta_post[:, None] * ct_grid[None, :]
    
            # pred: (n_reads, n_samples, n_ct)
            _pred = 1 / (1 + np.exp(-_kr_post[None, :, None] * (np.log(_reads) - _midpoint[None, :, :])))
            _pred_mean = _pred.mean(axis=1)[0]
            _pred_hdi = az.hdi(_pred[0].T, prob=0.95 )
    
            _ax.plot( ct_grid, _pred_mean, color='black', linewidth=1, linestyle='dashed', zorder=10 )
            _ax.fill_between( ct_grid, _pred_hdi[:,0], _pred_hdi[:,1], color='black', linewidth=0, alpha=0.1, zorder=5 )
    
            _points = _ax.scatter( "O1", "correct", data=vibe.loc[(vibe["sample_type"]==_st_name)&(vibe["reads"]==_reads)], zorder=20, alpha=0.4, linewidth=0)
            _points.set_clip_on(False)
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
            _ax.set_ylabel( f"Lineage w/\n{_mod} million reads", fontsize=10)
            _ax.grid( which="major", color="#efefef", linewidth=1 )
            _ax.grid( which="minor", color="#efefef", linewidth=0.5 )
            _ax.label_outer()
    
        plt.tight_layout()
        plt.savefig( f"analysis/plots/figureSX-lineage-vs-ct-{_mod}mil.pdf")
        plt.show()
    return


@app.cell
def _(calibration_curve, np, plt, vibe, vibes_trace):
    expected_p = vibes_trace.posterior_predictive["mu"].mean(dim=["chain", "draw"]).values
    observed_y = vibe["correct"].values

    # Calculate Brier Score (0 is perfect agreement, 0.25 is random guessing)
    brier_score = np.mean((expected_p - observed_y) ** 2)
    print(f"Brier Score: {brier_score:.4f}")

    prob_true, prob_pred = calibration_curve(observed_y, expected_p, n_bins=10)

    # Plotting the agreement
    plt.figure(figsize=(6, 6))
    plt.plot([0, 1], [0, 1], "k--", label="Perfect calibration")
    plt.plot(prob_pred, prob_true, "s-", label=f"Model (Brier: {brier_score:.3f})")
    plt.xlabel("Mean Predicted Probability")
    plt.ylabel("Fraction of True Successes")
    plt.title("Calibration Curve (Agreement Plot)")
    plt.legend()
    plt.show()
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
