from typing import *

import anndata as ad
import numpy as np
import pandas as pd

import plotly.graph_objects as go
from plotly.express.colors import sample_colorscale

import streamlit as st

from scrna_atlas.utils import collapse_polygons

cs = sample_colorscale("OrRd", list(np.linspace(0, 1, 50)))


def color_mapper(
    val: float,
    cs: List[str] = cs,
) -> str:
    """ """
    try:
        return cs[min(int(val * len(cs)), len(cs) - 1)]
    except:
        return cs[0]


def build_scatter(
    adata: ad.AnnData,
    color_umap_by: str,
    pick_umap_color: List[str],
    layer: Union[str, None] = "ntpm",
) -> List[go.Trace]:
    """ """
    scatter_data = []
    if "All" in pick_umap_color:
        for g in adata.obs[color_umap_by].unique():
            polygon_keys = adata.obs.loc[adata.obs[color_umap_by] == g].index.tolist()
            xpts, ypts = collapse_polygons(
                adata.uns["umap_polygons"][k] for k in polygon_keys
            )
            scatter_data.append(
                go.Scatter(
                    x=xpts,
                    y=ypts,
                    mode="lines",
                    marker={
                        "opacity": 0.25,
                    },
                    fill="toself",
                    name=g.upper(),
                )
            )
    else:
        # Categorical coloring of umap.
        if (
            len(pick_umap_color) == 0
            or pick_umap_color[0] in adata.obs[color_umap_by].tolist()
        ):
            # Pick points that should be colored gray.
            polygon_keys = adata.obs.loc[
                ~adata.obs[color_umap_by].isin(pick_umap_color)
            ].index.tolist()
            xpts, ypts = collapse_polygons(
                adata.uns["umap_polygons"][k] for k in polygon_keys
            )

            scatter_data.append(
                go.Scatter(
                    x=xpts,
                    y=ypts,
                    mode="lines",
                    marker={
                        "color": "grey",
                        "opacity": 0.25,
                    },
                    fill="toself",
                    name="Others",
                )
            )
            for g in adata.obs[color_umap_by].unique():
                if g in pick_umap_color:
                    polygon_keys = adata.obs.loc[
                        adata.obs[color_umap_by] == g
                    ].index.tolist()
                    xpts, ypts = collapse_polygons(
                        adata.uns["umap_polygons"][k] for k in polygon_keys
                    )
                    scatter_data.append(
                        go.Scatter(
                            x=xpts,
                            y=ypts,
                            mode="lines",
                            marker={
                                "opacity": 0.25,
                            },
                            fill="toself",
                            name=g.upper(),
                        )
                    )
        # Continuous coloring of umap.
        else:
            gs = adata.obs[color_umap_by].unique()
            obs = adata.obs.copy()
            for c in pick_umap_color:
                obs[c] = adata.layers[layer][:, adata.var["gene_name"] == c]
            exps = np.array(
                [
                    obs.loc[obs[color_umap_by] == g, pick_umap_color]
                    .mean(axis=0)
                    .values[0]
                    for g in gs
                ]
            )
            for g, e in zip(gs, exps):
                polygon_keys = adata.obs.loc[
                    adata.obs[color_umap_by] == g
                ].index.tolist()
                xpts, ypts = collapse_polygons(
                    adata.uns["umap_polygons"][k] for k in polygon_keys
                )
                scatter_data.append(
                    go.Scatter(
                        x=xpts,
                        y=ypts,
                        mode="lines",
                        line={
                            "width": 0.0,
                        },
                        marker={
                            "opacity": 1.0,
                            "color": color_mapper(e / exps.max()),
                        },
                        fill="toself",
                        name=g.upper(),
                    )
                )

    return scatter_data
