from pathlib import Path
from typing import List

import numpy as np
import pandas as pd
import anndata as ad
import streamlit as st
import plotly.graph_objects as go
from shapely.ops import unary_union
from shapely import Polygon, MultiPolygon

from scrna_atlas.multiadata import AtlasAdataWrapper

def collapse_polygons(
    exteriors: List[np.ndarray],
) -> List[np.ndarray]:
    """ """
    geoms = []
    for e in exteriors:
        geoms.append(Polygon(e.T))

    polygon = unary_union(geoms)
    
    x = []
    y = []

    if isinstance(polygon, MultiPolygon):
        for g in polygon.geoms:
            x.extend([p[0] for p in g.exterior.coords]+[None])
            y.extend([p[1] for p in g.exterior.coords]+[None])
    else:
        x.extend([p[0] for p in polygon.exterior.coords]+[None])
        y.extend([p[1] for p in polygon.exterior.coords]+[None])
        
    return x, y


def explorer_tab(
    multiadata: AtlasAdataWrapper,
    tab_title: str,
) -> None:
    st.title(tab_title)

    resolution_option = st.selectbox(
        'Pick atlas subset:',
        multiadata.adatas,
    )

    adata = multiadata.adatas[resolution_option]

    with st.sidebar:
        color_umap_by = st.selectbox(
            "Color umap by:",
            adata.obs.columns.unique(),
        )
        pick_umap_color = st.multiselect(
            "Color:",
            ['All']+adata.obs[color_umap_by].unique().tolist(),
        )

    scatter_data = []
    if 'All' in pick_umap_color:
        for g in adata.obs[color_umap_by].unique():
            polygon_keys = adata.obs.loc[adata.obs[color_umap_by] == g].index.tolist()
            xpts, ypts = collapse_polygons(
                adata.uns['umap_polygons'][k] for k in polygon_keys
            )
            scatter_data.append(
                # go.Scatter(
                #     x=adata.obsm['X_umap'][adata.obs[color_umap_by] == g,0][::30],
                #     y=adata.obsm['X_umap'][adata.obs[color_umap_by] == g,1][::30],
                #     mode='markers',
                #     marker={
                #         'opacity': 0.25,
                #     },
                #     name=g.upper(),
                # )
                
                go.Scatter(
                    # x=[p for k in polygon_keys for p in (adata.uns['umap_polygons'][k][0,:].tolist()+[None] if adata.uns['umap_polygons'][k].size > 0 else [None])],
                    # y=[p for k in polygon_keys for p in (adata.uns['umap_polygons'][k][1,:].tolist()+[None] if adata.uns['umap_polygons'][k].size > 0 else [None])],
                    x=xpts,
                    y=ypts,
                    mode='lines',
                    marker={
                        'opacity': 0.25,
                    },
                    fill='toself',
                    name=g.upper(),
                )
            )
    else:
        # Pick points that should be colored gray.
        scatter_data.append(
            go.Scatter(
                x=adata.obsm['X_umap'][~adata.obs[color_umap_by].isin(pick_umap_color),0][::30],
                y=adata.obsm['X_umap'][~adata.obs[color_umap_by].isin(pick_umap_color),1][::30],
                mode='markers',
                marker={
                    'color':'grey',
                    'opacity': 0.25,
                },
                name='Others',
            )
        )
        for g in adata.obs[color_umap_by].unique():
            if g in pick_umap_color:
                scatter_data.append(
                    go.Scatter(
                        x=adata.obsm['X_umap'][adata.obs[color_umap_by] == g,0][::30],
                        y=adata.obsm['X_umap'][adata.obs[color_umap_by] == g,1][::30],
                        mode='markers',
                        marker={
                            'opacity': 0.25,
                        },
                        name=g.upper(),
                    )
                )

    fig_umap = go.Figure(
        data=scatter_data,
    )
    fig_umap.update_layout(
        xaxis_title='umap_1',
        yaxis_title='umap_2',
        plot_bgcolor='rgba(0, 0, 0, 0)',
        paper_bgcolor='rgba(0, 0, 0, 0)',
    )

    st.plotly_chart(fig_umap, use_container_width=True, height=1500)
    with st.expander('Expression by Cell Type/Tissue'):
        celltype_tissue_select = st.radio("", ['Cell Type', 'Tissue'])
        gene_selection = st.multiselect(
            'Pick Gene/uProtein',
            adata.var.index,
            key=tab_title+'b'
        )
        bar_data = []

        for g in gene_selection:
            group_means = adata.obs.assign(
                gene_expression = adata.X[:,adata.var.index == g],
            ).groupby(
                'cell_type' if celltype_tissue_select == 'Cell Type' else 'tissue'
            ).agg('mean')
            st.write(group_means)
            # bar_data.append(
            #     go.Bar(
            #         name=g,
            #         x=adata.obs[
            #             'cell_type' if celltype_tissue_select == 'Cell Type' else 'tissue'
            #         ].unique(),
            #         y=
            #     )
            # )
            pass
        fig_expression = go.Figure(
            go.Bar(

            )
        )
        st.write(adata.obs['tissue'].unique())
    with st.expander('Genes'):
        st.write('Here')