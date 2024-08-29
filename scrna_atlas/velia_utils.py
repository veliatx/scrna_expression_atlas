import numpy as np
import pandas as pd
import streamlit as st

from pandas.api.types import (
    is_categorical_dtype,
    is_numeric_dtype,
)

from scrna_atlas import settings


@st.cache_data
def pull_sorfensembl_mappings(
    strip_versions: bool = True,
    filter_secreted: bool = True,
    filter_riboseq: bool = True,
) -> None:
    """ """
    sorf_df = pd.read_parquet(settings.DASHBOARD_SORFS_URL)

    if filter_secreted:
        sorf_df = filter_sorfensembl_mapping_secreted(sorf_df)
    if filter_riboseq:
        sorf_df = filter_sorfensembl_mapping_riboseq(sorf_df)
        sorf_df = sorf_df.loc[sorf_df["Ribo-Seq sORF"]]

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


def filter_sorfensembl_mapping_secreted(
    df: pd.DataFrame,
) -> pd.DataFrame:
    """ """
    signal_cols = [
        "SignalP 4.1_cut",
        "SignalP 5b_cut",
        "SignalP 6slow_cut",
        "Deepsig_cut",
    ]
    measured_secreted_or_predicted_secreted = df["secreted"] | (
        df[signal_cols] > -1
    ).any(axis=1)

    return df.loc[measured_secreted_or_predicted_secreted]


def filter_sorfensembl_mapping_riboseq(df):
    """
    Temporary function to enforce ribo-seq filtering
    for primary entries in collection
    """
    df["Ribo-Seq sORF"] = (
        (df["source"].apply(lambda x: "gencode_riboseq" in x))
        | (df["source"].apply(lambda x: "velia_phase1_Bona fide" in x))
        | (df["source"].apply(lambda x: "velia_phase1_Chen" in x))
        | (df["source"].apply(lambda x: "velia_phase1_Prensner" in x))
        | (df["source"].apply(lambda x: "velia_phase2_Chang_Saghatelian" in x))
        | (df["source"].apply(lambda x: "velia_phase2_Chothani2022_SignalP" in x))
        | (df["source"].apply(lambda x: "velia_phase2_Bianca_Chen" in x))
        | (df["source"].apply(lambda x: "velia_phase2_Bonafide_Bianca" in x))
        | (df["source"].apply(lambda x: "velia_phase2_Cao_Slavoff_MINAS60" in x))
        | (df["source"].apply(lambda x: "velia_phase2_Rat_Cardiac_Huang" in x))
        | (df["source"].apply(lambda x: "velia_phase2_Mudge2022_SignalP" in x))
        | (df["source"].apply(lambda x: "velia_phase5_Blume_Mudge" in x))
        | (df["source"].apply(lambda x: "velia_phase5_bona fide" in x))
        | (df["source"].apply(lambda x: "velia_phase6_plasma_mass_spec" in x))
        | (df["source"].apply(lambda x: "velia_phase6_public_mass_spec" in x))
        | (df["source"].apply(lambda x: "velia_phase9_orfrater" in x))
        | (df["source"].apply(lambda x: "velia_phase9_Olsen" in x))
        | (df["source"].apply(lambda x: "velia_phase9_Li et al VSMC" in x))
        | (df["source"].apply(lambda x: "velia_phase7_Ribo-seq_PBMC_LPS_R848" in x))
        | (df["source"].apply(lambda x: "velia_phase10_riboseq_230114" in x))
        | (df["source"].apply(lambda x: "velia_phase11_riboseq_240214" in x))
        | (df["source"].apply(lambda x: "swissprot" in x))
        | (df["source"].apply(lambda x: "MetaORF v1.0" in x))
        | (df["source"].apply(lambda x: "ENSEMBL" in x))
        | (df["source"].apply(lambda x: "BestRefSeq" in x))
        | (df["source"].apply(lambda x: "velia_phase5_uniprot-tremble" in x))
        | (df["screening_phase"] == "Not Screened")
        | (df["orf_xrefs"].astype(str).str.contains("RibORF"))
        | (df["protein_xrefs"].astype(str).str.contains("RibORF"))
    )

    ribo_df = df[df["Ribo-Seq sORF"]].copy()
    x = ribo_df.groupby("aa").agg(list)

    vtx_to_keep = []

    for i, row in x.iterrows():
        vtx_id = ""

        if len(row["vtx_id"]) > 1:

            for j, phase in enumerate(row["screening_phase"]):
                if "phase" in phase.lower():
                    vtx_id = row["vtx_id"][j]

            if vtx_id == "":
                vtx_id = row["vtx_id"][0]
        else:
            vtx_id = row["vtx_id"][0]

        vtx_to_keep.append(vtx_id)

    ribo_df = ribo_df[ribo_df["vtx_id"].isin(vtx_to_keep)].copy()

    ribo_aa = set(ribo_df["aa"])

    non_ribo_df = df[~df["Ribo-Seq sORF"]].copy()
    non_ribo_df = non_ribo_df[~non_ribo_df["aa"].isin(ribo_aa)]

    df = pd.concat([ribo_df, non_ribo_df])

    return df


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
        .unique()
    )
    return pd.DataFrame(
        np.ones(ensembl_id.size, dtype=bool),
        columns=["receptor"],
        index=ensembl_id,
    )


def filter_dataframe_dynamic(df: pd.DataFrame) -> pd.DataFrame:
    """
    Adds a UI on top of a dataframe to let viewers filter columns

    Args:
        df (pd.DataFrame): Original dataframe

    Returns:
        pd.DataFrame: Filtered dataframe
    """

    df = df.copy()

    modification_container = st.container()

    with modification_container:
        to_filter_columns = st.multiselect("Filter dataframe on", df.columns)
        for column in to_filter_columns:
            left, right = st.columns((1, 20))
            if is_numeric_dtype(df[column]):
                _min = float(df[column].min())
                _max = float(df[column].max())
                step = (_max - _min) / 100
                user_num_input = right.slider(
                    f"Values for {column}",
                    min_value=_min,
                    max_value=_max,
                    value=(_min, _max),
                    step=step,
                )
                df = df[df[column].between(*user_num_input)]
            else:
                user_text_input = right.text_input(
                    f"Substring or regex in {column}",
                )
                if user_text_input:
                    df = df[df[column].astype(str).str.contains(user_text_input)]

    return df
