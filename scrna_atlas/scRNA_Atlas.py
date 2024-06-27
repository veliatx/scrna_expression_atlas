from pathlib import Path

import plotly.graph_objects as go
import pandas as pd
import anndata as ad
import streamlit as st


def main():
    st.set_page_config(layout="wide")

    st.title("scRNA atlas")
    st.write("Pick existing scRNA atlas or individual experiments.")
    st.write(
        """Data come from 3-4 sources:
1. Tabula Sapiens -> scRNA-seq from CellXGene tabula sapiens.
2. HPA scRNA -> scRNA-seq from Human Protein Atlas. 
3. Tabula Muris -> scRNA-seq from CellXGene tabula muris.
4. Experiment Explorer -> Individual CellXGene experiments or external scRNA-seq reprocessed to CellXGene standards.
"""
    )


if __name__ == "__main__":
    main()
