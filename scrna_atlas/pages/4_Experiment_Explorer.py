import streamlit as st

from scrna_atlas.atlas_explorer import explorer_tab

def main():
    select_experiment = st.selectbox(
        'Select experiment to explor:',
        ['A', 'B', 'C', 'D'],
    )

if __name__ == '__main__':
    main()