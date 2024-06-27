from pathlib import Path

import streamlit as st

from scrna_atlas.atlas_explorer import atlas_explorer_tab
from scrna_atlas.multiadata import AtlasAdataWrapper


def main():
    st.set_page_config(layout="wide")

    if "multiadata_muris" not in st.session_state:
        multiadata = AtlasAdataWrapper(
            Path("data/tabula_muris/10X"),
            "tabula_muris",
            full_atlas="All - A single-cell transcriptomic atlas characterizes ageing tissues in the mouse - 10x",
        )
        st.session_state["multiadata_muris"] = multiadata
    else:
        multiadata = st.session_state["multiadata_muris"]

    atlas_explorer_tab(
        multiadata,
        "Tabula Muris",
    )


if __name__ == "__main__":
    main()
