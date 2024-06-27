from pathlib import Path

import streamlit as st

from scrna_atlas.atlas_explorer import atlas_explorer_tab
from scrna_atlas.multiadata import AtlasAdataWrapper


def main():
    st.set_page_config(layout="wide")

    if "multiadata_sapiens" not in st.session_state:
        multiadata = AtlasAdataWrapper(
            Path("data/tabula_sapiens/10X"),
            "tabula_sapiens",
            full_atlas="Tabula Sapiens - All Cells",
        )
        st.session_state["multiadata_sapiens"] = multiadata
    else:
        multiadata = st.session_state["multiadata_sapiens"]
    atlas_explorer_tab(
        multiadata,
        "Tabula Sapiens",
    )


if __name__ == "__main__":
    main()
