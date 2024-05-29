from pathlib import Path 

import streamlit as st

from scrna_atlas.atlas_explorer import explorer_tab
from scrna_atlas.multiadata import AtlasAdataWrapper

def main():
    multiadata = AtlasAdataWrapper(
        Path('data/tabula_sapiens'),
        'tabula_sapiens',
    )
    explorer_tab(
        multiadata,
        'Tabula Sapiens',
    )


if __name__ == '__main__':
    main()