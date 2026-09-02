import marimo

__generated_with = "0.23.16"
app = marimo.App()


@app.cell
def _():
    import marimo as mo
    import pandas as pd
    import numpy as np

    samples = pd.read_csv( "resources/samples.csv")
    return np, pd, samples


@app.cell
def _(samples):
    samples.head()
    return


@app.cell
def _(np, pd):
    oct = pd.read_csv( "resources/WGS-Sample-type-QPCR.csv" )
    oct = oct.rename( columns={"Sample ID" : "sample", "WS Qiagen" : "sample_type"} )
    oct["sample_type"] = oct["sample_type"].str.slice(5).replace({
        "RDTWS_APW" : "RDT_WS_APW",
        "FPWS_APW" : "FP_WS_APW",
        "RDTIS_TCBS" : "RDT_IS_TCBS",
        "FPIS_TCBS" : "FP_IS_TCBS",
    })
    oct["toxR"] = oct["toxR"].replace( {"XX": np.nan, "0" : np.nan} ).astype(float)
    oct["ctxA"] = oct["ctxA"].replace( {"XX": np.nan, "0" : np.nan} ).astype(float)
    oct
    return (oct,)


@app.cell
def _(oct, pd, samples):
    ct = pd.read_csv( "resources/WGS QPCR Combined Results_O1.csv")
    ct.columns = ["sample", "sample_type", "target", "O1"]
    ct = ct.drop( columns=["target"])
    ct["sample"] = ct["sample"].str.strip()
    ct["O1"] = ct["O1"].replace( {"Undetermined" : "40"} ).astype( float )
    ct["sample_type"] = ct["sample_type"].str.replace( " ", "_" ).replace( {
        "WSFP" : "WS_FP",
        "WSQ" : "WS_Q",
        "WSRDT" : "WS_RDT",
    }).str.strip()
    ct = ct.merge( oct, on=["sample", "sample_type"], how="outer" )
    ct = ct.merge( samples[["cell", "sample", "sample_type"]], on=["sample", "sample_type"], how="left" )
    ct = ct.reindex(columns=["cell", "sample", "sample_type", "O1", "toxR", "ctxA"])
    ct
    return (ct,)


@app.cell
def _(ct):
    ct.to_csv( "resources/WGS-qPCR.csv", index=False )
    return


if __name__ == "__main__":
    app.run()
