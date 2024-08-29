from pathlib import Path
from datetime import datetime
from typing import *

import numpy as np
import pandas as pd
import anndata as ad
import streamlit as st

from plotly.subplots import make_subplots
import plotly.graph_objects as go

from scrna_atlas import settings, velia_utils
from scrna_atlas.plotting import build_scatter


def experiment_explorer_tab(
    adata: ad.AnnData,
    tab_title: str,
    full_adata: Union[ad.AnnData, None],
) -> None:
    """ """
    st.header(tab_title)
    meta = adata.uns["meta"]
    st.markdown(
        """##### Experiment metadata:
* Species: {species} 
* Pubmed ID: {pmid}
* Assay: {assay}
* Comments: {comments}
""".format(
            species=meta.get("species", None),
            pmid=meta.get("pmid", None),
            assay=",".join(meta.get("assay", "")),
            comments=meta.get("comments", None),
        )
    )

    with st.sidebar:
        columns_display = adata.obs.columns[
            adata.obs.columns.isin(
                [
                    "cell_type",
                    "tissue",
                    "donor_id",
                    "tissue_cell_type",
                    "assay",
                    "disease",
                    "louvain",
                    "leiden",
                    "sex",
                    "cluster",
                    "cell_type_cluster",
                ]
                + adata.uns["factors"].tolist()
            )
        ]
        color_umap_by = st.selectbox(
            "Color umap by:",
            columns_display,
            index=(
                columns_display.get_loc("cell_type")
                if "cell_type" in columns_display
                else 0
            ),
        )
        pick_umap_color = st.multiselect(
            "Color:",
            ["All"] + adata.obs[color_umap_by].unique().tolist(),
            default=["All"],
        )
    if full_adata:
        umap_plot_polygons = st.toggle("Polygons", value=True, key="umap_plot_polygons")
    else:
        umap_plot_polygons = True
    scatter_data = build_scatter(adata, color_umap_by, pick_umap_color)
    fig_umap = go.Figure(data=scatter_data)
    fig_umap.update_layout(
        xaxis_title="umap_1",
        yaxis_title="umap_2",
        plot_bgcolor="rgba(0, 0, 0, 0)",
        paper_bgcolor="rgba(0, 0, 0, 0)",
    )

    st.plotly_chart(fig_umap, use_container_width=True, height=1500)

    with st.expander("UMAPs"):
        if full_adata:
            umap_plot_polygons = st.toggle(
                "Polygons", value=True, key="umap_compare_plot_polygons"
            )
        else:
            umap_plot_polygons = True

        columns_display = adata.obs.columns[
            adata.obs.columns.isin(
                [
                    "cell_type",
                    "tissue",
                    "donor_id",
                    "tissue_cell_type",
                    "assay",
                    "disease",
                    "louvain",
                    "leiden",
                    "sex",
                    "cluster",
                    "cell_type_cluster",
                ]
                + adata.uns["factors"].tolist()
            )
        ]
        umapexp_color_umap_by = st.multiselect(
            "Color umap by (up to 3):",
            columns_display,
        )
        if len(umapexp_color_umap_by) > 3:
            st.write("Select no more than 3 metadata.")
        elif len(umapexp_color_umap_by) > 0:
            fig_umapexp = make_subplots(
                rows=1,
                cols=len(umapexp_color_umap_by),
                shared_xaxes=True,
                shared_yaxes=True,
                subplot_titles=umapexp_color_umap_by,
            )
            for i, c in enumerate(umapexp_color_umap_by):
                c_scatter_data = build_scatter(adata, c, ["All"])
                for cs in c_scatter_data:
                    fig_umapexp.add_trace(cs, row=1, col=i + 1)
                    fig_umapexp.update_layout({f"xaxis{i+1}_title": "umap_1"})
            fig_umapexp.update_layout(
                xaxis_title="umap_1",
                yaxis_title="umap_2",
                plot_bgcolor="rgba(0, 0, 0, 0)",
                paper_bgcolor="rgba(0, 0, 0, 0)",
                showlegend=False,
            )
            st.plotly_chart(fig_umapexp, use_container_width=True, height=1000)

    with st.expander("Expression"):
        celltype_tissue_select = st.multiselect(
            "Expression by Cell Type/Factor:",
            ["Cell Type"]
            + list(map(lambda x: x.capitalize(), adata.uns["factors"]))
            + list(
                map(
                    lambda x: x.capitalize(),
                    adata.obs.columns[
                        adata.obs.columns.isin(
                            ["sex", "tissue", "leiden", "cell_type_cluster"]
                        )
                    ],
                )
            ),
        )
        celltype_tissue_select = list(
            map(lambda x: x.lower().replace(" ", "_"), celltype_tissue_select)
        )

        var_wvtx = adata.var.merge(
            velia_utils.pull_sorfensembl_mappings(),
            left_index=True,
            right_index=True,
            how="left",
        ).assign(
            vtx_id_concat=lambda x: x["vtx_id"].map(
                lambda y: ",".join(y) if isinstance(y, list) else "NA"
            )
        )
        gene_selection_options = np.concatenate(
            [
                var_wvtx["gene_name"].values,
                var_wvtx["vtx_id"].explode().drop_duplicates().values,
            ]
        )

        gene_selection = st.selectbox(
            "Pick Gene/uProtein",
            gene_selection_options,
            key=tab_title + "b",
        )
        log_yaxis_select = st.radio("Scale", ["linear", "log"], horizontal=True)

        if len(gene_selection) > 0 and len(celltype_tissue_select) > 0:
            gene_name_locs = np.where(var_wvtx["gene_name"] == gene_selection)[0]
            vtx_name_locs = np.where(
                var_wvtx["vtx_id_concat"].map(lambda x: gene_selection in x)
            )[0]
            column_names = (
                var_wvtx["gene_name"][gene_name_locs.tolist()].tolist()
                + var_wvtx["vtx_id_concat"][vtx_name_locs.tolist()].tolist()
            )
            gene_selection_df = pd.concat(
                [
                    adata.obs.copy(),
                    pd.DataFrame(
                        adata.layers[
                            "ntpm" if "ntpm" in adata.layers else "count_mean"
                        ][:, gene_name_locs.tolist() + vtx_name_locs.tolist()],
                        index=adata.obs.index,
                        columns=column_names,
                    ),
                ],
                axis=1,
                ignore_index=False,
            )

            box_data = go.Box(
                name=gene_selection,
                x=[
                    "__".join(r)
                    for i, r in gene_selection_df.loc[
                        :, celltype_tissue_select
                    ].iterrows()
                ],
                y=gene_selection_df[column_names[0]].tolist(),
                boxpoints="all",
                fillcolor="grey",
                line_color="darkgrey",
            )

            fig_expression = go.Figure(
                data=box_data,
            )
            fig_expression.update_layout(
                plot_bgcolor="rgba(0, 0, 0, 0)",
                paper_bgcolor="rgba(0, 0, 0, 0)",
                yaxis_type="log" if log_yaxis_select == "log" else "linear",
                yaxis_title="log nTPM" if log_yaxis_select == "log" else "nTPM",
            )
            st.plotly_chart(fig_expression, use_container_width=True, height=1500)

    with st.expander("Differential Expression"):
        if not any([k.startswith("de-") for k in adata.varm.keys()]):
            st.write("No differential expression available for this experiment.")
        else:
            contrast_selection = st.selectbox(
                "Pick contrast:",
                [k for k in adata.varm.keys() if k.startswith("de-")],
                format_func=lambda x: x.split("de-")[-1],
            )
            # Fill na with 1.0, na likely because of independent filtering step of de pipeline.
            contrast_selection_df = (
                adata.varm[contrast_selection]
                .loc[~adata.varm[contrast_selection]["log2FoldChange"].isna()]
                .fillna({"pvalue": 1.0, "padj": 1.0})
                .merge(
                    velia_utils.pull_sorfensembl_mappings(),
                    left_index=True,
                    right_index=True,
                    how="left",
                )
                .merge(
                    velia_utils.pull_receptorensembl_mappings(),
                    left_index=True,
                    right_index=True,
                    how="left",
                )
                .assign(
                    vtx_id=lambda x: x.vtx_id.map(
                        lambda y: y if not isinstance(y, float) else []
                    )
                )
            )

            st.markdown(
                """###### Total genes tested: {n}\n###### Pseudobulk info:""".format(
                    n=contrast_selection_df.shape[0]
                )
            )
            if "tissue_cell_type" in adata.obs.columns:
                obs_filter = (
                    adata.obs["cell_type"].str.replace(" ", "")
                    == contrast_selection.split("__")[1].split("de-")[-1]
                )
            elif "cell_type_cluster" in adata.obs.columns:
                obs_filter = (
                    adata.obs["cell_type_cluster"].str.replace(" ", "")
                    == contrast_selection.split("__")[0].split("de-")[-1]
                )
            else:
                obs_filter = (
                    adata.obs["cell_type"].str.replace(" ", "")
                    == contrast_selection.split("__")[0].split("de-")[-1]
                )
            st.dataframe(
                adata.obs.loc[
                    obs_filter,
                    adata.obs.columns[
                        adata.obs.columns.isin(
                            [
                                "bulk_n_cells",
                                "bulk_total_cells",
                                "cell_type",
                                "tissue",
                                "donor_id",
                                "sex",
                            ]
                            + adata.uns["factors"].tolist()
                        )
                    ],
                ].rename({"bulk_total_cells": "n_reads"}, axis=1)
            )
            de_col1, de_col2 = st.columns(2)
            with de_col1:
                st.subheader("Plots:")
                volcano_tab, ma_tab = st.tabs(["Volcano:", "MA:"])
                with volcano_tab:
                    plot_fig = go.Figure(
                        go.Scatter(
                            x=contrast_selection_df["log2FoldChange"],
                            y=contrast_selection_df["log10_padj"],
                            mode="markers",
                            marker={
                                "opacity": 0.25,
                                "color": list(
                                    map(
                                        lambda x: "blue" if x else "red",
                                        (
                                            np.abs(
                                                contrast_selection_df["log2FoldChange"]
                                            )
                                            > settings.DE_LFC_THRESH
                                        )
                                        & (
                                            contrast_selection_df["padj"]
                                            < settings.DE_PADJ_THRESH
                                        ),
                                    )
                                ),
                            },
                            customdata=contrast_selection_df[
                                ["gene_name", "log2FoldChange", "padj", "vtx_id"]
                            ],
                            hovertemplate="<b>Gene: %{customdata[0]}</b><br>"
                            + "<b>log2FC: %{customdata[1]:,.4f}</b><br>"
                            + "<b>padj: %{customdata[2]:,.4f}</b><br>"
                            + "<b>uPs: %{customdata[3]}</b>",
                        ),
                    )
                    plot_fig.update_layout(
                        yaxis_autorange="reversed",
                        xaxis_title="log2FoldChange",
                        yaxis_title="log10_padj",
                        plot_bgcolor="rgba(0, 0, 0, 0)",
                        paper_bgcolor="rgba(0, 0, 0, 0)",
                    )
                    st.plotly_chart(plot_fig, use_container_width=True, height=1500)
                with ma_tab:
                    plot_fig = go.Figure(
                        go.Scatter(
                            x=np.log2(contrast_selection_df["baseMean"]),
                            y=contrast_selection_df["log2FoldChange"],
                            mode="markers",
                            marker={
                                "opacity": 0.25,
                                "color": list(
                                    map(
                                        lambda x: "blue" if x else "red",
                                        (
                                            np.abs(
                                                contrast_selection_df["log2FoldChange"]
                                            )
                                            > settings.DE_LFC_THRESH
                                        )
                                        & (
                                            contrast_selection_df["padj"]
                                            < settings.DE_PADJ_THRESH
                                        ),
                                    )
                                ),
                            },
                            customdata=contrast_selection_df[
                                ["gene_name", "log2FoldChange", "padj", "vtx_id"]
                            ],
                            hovertemplate="<b>Gene: %{customdata[0]}</b><br>"
                            + "<b>log2FC: %{customdata[1]:,.4f}</b><br>"
                            + "<b>padj: %{customdata[2]:,.4f}</b><br>"
                            + "<b>uPs: %{customdata[3]}</b>",
                        ),
                    )
                    plot_fig.update_layout(
                        xaxis_title="log2AveExpression",
                        yaxis_title="log2FoldChange",
                        plot_bgcolor="rgba(0, 0, 0, 0)",
                        paper_bgcolor="rgba(0, 0, 0, 0)",
                    )
                    st.plotly_chart(plot_fig, use_container_width=True, height=1500)
            with de_col2:
                st.subheader("DE List:")
                st.dataframe(
                    contrast_selection_df.loc[
                        :,
                        ~contrast_selection_df.columns.isin(
                            ["baseMean", "lfcSE", "stat", "log10_padj"]
                        ),
                    ].sort_values("padj")
                )
                st.download_button(
                    label="Download DE List",
                    data=contrast_selection_df.to_csv(),
                    file_name=f'{tab_title.replace(" ","")}_{contrast_selection}_{datetime.now().strftime("%Y%m%d")}.csv',
                    mime="text/csv",
                )
            var_wvtx = adata.var.merge(
                velia_utils.pull_sorfensembl_mappings(),
                left_index=True,
                right_index=True,
                how="left",
            ).assign(
                vtx_id_concat=lambda x: x["vtx_id"].map(
                    lambda y: ",".join(y) if isinstance(y, list) else "NA"
                )
            )
            gene_selection_options = np.concatenate(
                [
                    var_wvtx["gene_name"].values,
                    var_wvtx["vtx_id"].explode().drop_duplicates().values,
                ]
            )
            st.subheader("Seach contrasts by gene/uP:")
            search_gene = st.selectbox(
                "Search contrasts by gene/uP",
                gene_selection_options,
                label_visibility="collapsed",
            )

            search_gene_mask = (var_wvtx["gene_name"] == search_gene) | var_wvtx[
                "vtx_id_concat"
            ].map(lambda x: search_gene in x)

            search_gene_button = st.button("Search")
            if search_gene_button:
                search_results = []
                for k, v in adata.varm.items():
                    v_series = v.loc[search_gene_mask]
                    if v_series.size == 0:
                        continue
                    v_series["contrast"] = k
                    if np.isfinite(v_series["log2FoldChange"].values[0]):
                        search_results.append(v_series)
                if len(search_results) > 0:
                    search_df = (
                        pd.concat(search_results)
                        .merge(
                            velia_utils.pull_sorfensembl_mappings(),
                            left_index=True,
                            right_index=True,
                            how="left",
                        )
                        .merge(
                            velia_utils.pull_receptorensembl_mappings(),
                            left_index=True,
                            right_index=True,
                            how="left",
                        )
                    )
                    st.dataframe(
                        search_df.loc[
                            :,
                            ~search_df.columns.isin(
                                ["baseMean", "lfcSE", "stat", "log10_padj"]
                            ),
                        ]
                    )
                else:
                    st.write(
                        "Gene not expressed/Not enough samples to perform DE testing..."
                    )
