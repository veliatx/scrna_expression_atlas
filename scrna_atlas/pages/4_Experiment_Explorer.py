import pandas as pd
import streamlit as st
import anndata as ad

from scrna_atlas.experiment_explorer import experiment_explorer_tab
from scrna_atlas import settings


def main():
    st.set_page_config(layout="wide")

    if "experiment_df" not in st.session_state:
        experiment_df = pd.read_csv(settings.EXPERIMENT_EXPLORER_INDEX, sep="\t")
        st.session_state["experiment_df"] = experiment_df.copy()
    else:
        experiment_df = st.session_state["experiment_df"]

    select_experiment = st.selectbox(
        "Select experiment to explore:",
        experiment_df["title"].tolist(),
    )

    if (
        "active_de_adata" not in st.session_state.keys()
        or st.session_state["active_de_adata"][0] != select_experiment
    ):
        active_padata = ad.read_h5ad(
            settings.EXPERIMENT_LOC
            / experiment_df.loc[
                experiment_df["title"] == select_experiment, "pseudobulk_path"
            ].values[0],
            backed="r",
        )
        adata_name = experiment_df.loc[
            experiment_df["title"] == select_experiment, "full_path"
        ].values[0]
        adata_path = (
            (settings.EXPERIMENT_LOC / adata_name)
            if adata_name != "None" and not isinstance(adata_name, float)
            else None
        )
        active_adata = (
            adata_path if not adata_path else ad.read_h5ad(adata_path, backed="r")
        )
        st.session_state["active_de_adata"] = (
            select_experiment,
            active_padata,
            active_adata,
        )
    else:
        _, active_padata, active_adata = st.session_state["active_de_adata"]

    experiment_explorer_tab(
        active_padata,
        select_experiment,
        full_adata=active_adata,
    )


if __name__ == "__main__":
    main()
