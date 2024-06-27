import numpy as np
import pandas as pd
import streamlit as st

from scrna_atlas import settings


@st.cache_data
def pull_sorfensembl_mappings(
    strip_versions: bool = True,
) -> None:
    """ """
    sorf_df = pd.read_parquet(settings.DASHBOARD_SORFS_URL)
    sorf_df["gene_xrefs"] = sorf_df["gene_xrefs"].str.split(";")
    sorf_df = (
        sorf_df.reset_index(drop=True)
        .loc[:, ["vtx_id", "gene_xrefs"]]
        .explode("gene_xrefs")
    )
    if strip_versions:
        sorf_df["gene_xrefs"] = sorf_df["gene_xrefs"].str.split(".").str[0]
    sorf_df = (
        sorf_df.loc[sorf_df["gene_xrefs"].str.startswith("ENSG")]
        .groupby("gene_xrefs")[["vtx_id"]]
        .agg(list)
    )

    return sorf_df


@st.cache_data
def pull_receptorensembl_mappings():
    """Parsing ensembl ids out of HPA right now. Better to map between ensemlb transcript
    id to ensembl gene_id, but this works for now.
    """
    receptor_df = pd.read_csv(settings.RECEPTORS_URL, converters={"xrefs": pd.eval})
    receptor_df = receptor_df.loc[
        receptor_df["Receptor"]
        & receptor_df["xrefs"].map(lambda x: any(["HPA:ENSG" in xr for xr in x]))
    ]
    ensembl_id = (
        receptor_df["xrefs"]
        .map(lambda x: [xr.split("HPA:")[-1] for xr in x if xr.startswith("HPA:")][0])
        .values
    )
    return pd.DataFrame(
        np.ones(ensembl_id.size, dtype=bool),
        columns=["receptor"],
        index=ensembl_id,
    )
