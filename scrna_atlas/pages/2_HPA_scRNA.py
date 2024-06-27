from pathlib import Path

import streamlit as st

from scrna_atlas.atlas_explorer import atlas_explorer_tab
from scrna_atlas.multiadata import AtlasAdataWrapper


def main():
    st.set_page_config(layout="wide")

    if "multiadata_hpa" not in st.session_state:
        multiadata = AtlasAdataWrapper(
            Path("data/hpa/"),
            "hpa",
            glob_pattern="*/*.h5ad",
            full_atlas="scRNA HPA - All",
        )
        st.session_state["multiadata_hpa"] = multiadata
    else:
        multiadata = st.session_state["multiadata_hpa"]

    atlas_explorer_tab(
        multiadata,
        "HPA scRNA",
    )


if __name__ == "__main__":
    main()
