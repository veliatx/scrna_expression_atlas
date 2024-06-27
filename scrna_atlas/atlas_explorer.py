from pathlib import Path
from typing import List
import warnings

import numpy as np
import pandas as pd
import anndata as ad
import streamlit as st

import plotly.graph_objects as go
from plotly.subplots import make_subplots

from shapely.ops import unary_union
from shapely import Polygon, MultiPolygon, GeometryCollection

from scrna_atlas.multiadata import AtlasAdataWrapper
from scrna_atlas import velia_utils, utils
from scrna_atlas.plotting import build_scatter


def atlas_explorer_tab(
    multiadata: AtlasAdataWrapper,
    tab_title: str,
) -> None:
    st.title(tab_title)

    resolution_option = st.selectbox(
        "Pick atlas subset:",
        multiadata.list_adatas,
    )

    adata = multiadata.adatas[resolution_option]

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
                    "sex",
                    "cluster",
                ]
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

    meta = adata.uns.get("meta", None)
    if not meta:
        meta = {}
    st.markdown(
        """##### Comments:
"""
        + " ".join(meta.get("comments", []))
    )

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
                    "sex",
                    "cluster",
                ]
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
        celltype_tissue_select = st.radio(
            "Expression by Cell Type/Tissue",
            ["Cell Type", "Tissue", "Tissue + Cell Type", "Louvain"],
            horizontal=True,
        )
        if celltype_tissue_select == "Tissue + Cell Type":
            celltype_tissue_select = "Tissue Cell Type"
        celltype_tissue_select = celltype_tissue_select.lower().replace(" ", "_")

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
        gene_selection = st.multiselect(
            "Pick Gene/uProtein", gene_selection_options, key=tab_title + "b"
        )
        log_yaxis_select = st.radio("Scale", ["linear", "log"], horizontal=True)
        bar_data = []

        gene_name_locs = np.where(var_wvtx["gene_name"].isin(gene_selection))[0]
        vtx_name_locs = np.where(
            var_wvtx["vtx_id_concat"].map(lambda x: any(g in x for g in gene_selection))
        )[0]
        column_names = (
            var_wvtx["gene_name"][gene_name_locs.tolist()].tolist()
            + var_wvtx["vtx_id_concat"][vtx_name_locs.tolist()].tolist()
        )
        gene_selection_df = pd.concat(
            [
                adata.obs.copy(),
                pd.DataFrame(
                    adata.layers["ntpm"][
                        :, gene_name_locs.tolist() + vtx_name_locs.tolist()
                    ],
                    index=adata.obs.index,
                    columns=column_names,
                ),
            ],
            axis=1,
            ignore_index=False,
        )

        for g in column_names:
            g_ix = [i for i, s in enumerate(gene_selection) if s in g][0]
            gene_selection_df_groupby = gene_selection_df.groupby(
                celltype_tissue_select
            )
            group_means = gene_selection_df_groupby[g].agg(np.mean).fillna(0.0)
            group_ncells = (
                gene_selection_df_groupby["bulk_n_cells"].agg(np.sum).fillna(0.0)
            )
            bar_data.append(
                go.Bar(
                    name=gene_selection[g_ix],
                    x=group_means.index,
                    y=group_means,
                    text=group_ncells,
                )
            )

        fig_expression = go.Figure(
            data=bar_data,
        )
        fig_expression.update_layout(
            plot_bgcolor="rgba(0, 0, 0, 0)",
            paper_bgcolor="rgba(0, 0, 0, 0)",
            yaxis_title="log nTPM" if log_yaxis_select == "log" else "nTPM",
            yaxis_type="log" if log_yaxis_select == "log" else "linear",
            xaxis_categoryorder="total descending",
            xaxis_tickmode="linear",
        )
        st.plotly_chart(fig_expression, use_container_width=True, height=1500)
    with st.expander("Cell Type/Tissue Specificity"):
        tau_adata = multiadata.full_atlas_adata

        specificity_options = ["Cell Type", "Tissue", "Tissue + Cell Type", "Louvain"]
        if "cluster" in tau_adata.obs.columns:
            specificity_options.append("Cluster")
        specificity_by = st.radio(
            "Filter specificity by:",
            specificity_options,
            horizontal=True,
        )
        if specificity_by == "Tissue + Cell Type":
            specificity_by = "Tissue Cell Type"
        specificity_by = specificity_by.lower().replace(" ", "_")

        # tau_or_label_specificity = st.toggle("Tau/Label Specificity:", value=False)
        tau_or_label_specificity = False
        if (
            not tau_or_label_specificity
            and tau_adata.var.columns.str.startswith("tau").any()
        ):
            tau_df = (
                tau_adata.var.loc[
                    :,
                    (
                        ~tau_adata.var.columns.str.startswith("tau")
                        | (tau_adata.var.columns == f"tau__{specificity_by}")
                    )
                    & (
                        ~tau_adata.var.columns.isin(
                            [
                                "n_cells",
                                "means",
                                "dispersions",
                                "dispersions_norm",
                                "highly_variable",
                                "feature_is_filtered",
                                "feature_name",
                                "feature_reference",
                                "feature_biotype",
                                "feature_length",
                                "feature_type",
                                "mean",
                                "std",
                            ]
                        )
                    ),
                ]
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
                .rename({f"tau__{specificity_by}": "tau"}, axis=1)
            )
            st_col1, st_col2 = st.columns(2)
            with st_col1:
                fig_tau = go.Figure(
                    go.Histogram(
                        x=tau_df["tau"],
                        nbinsx=50,
                    )
                )
                fig_tau.update_layout(
                    plot_bgcolor="rgba(0, 0, 0, 0)",
                    paper_bgcolor="rgba(0, 0, 0, 0)",
                    xaxis_title="tau",
                )

                st.plotly_chart(fig_tau, use_container_width=True, height=1000)

                tau_select = st.slider(
                    "Select tau:",
                    0.0,
                    1.0,
                    (0.7, 1.0),
                )
            with st_col2:
                display_tau_df = tau_df.loc[
                    (tau_df["tau"] >= tau_select[0])
                    & (tau_df["tau"] <= tau_select[1])
                    & ~tau_df["tau"].isna(),
                ].sort_values("tau", ascending=False)
                st.write(
                    f"Genes: {display_tau_df.shape[0]}, uPs: {display_tau_df['vtx_id'].explode().shape[0]}"
                )
                tau_select = st.dataframe(
                    display_tau_df,
                    use_container_width=True,
                    selection_mode="single-row",
                    on_select="rerun",
                )["selection"]["rows"]
            if tau_select:
                st_select_col1, st_select_col2 = st.columns(2)
                with st_select_col1:
                    fig_tau_scatter1 = go.Figure(
                        data=build_scatter(
                            tau_adata,
                            color_umap_by=specificity_by,
                            pick_umap_color=["All"],
                        )
                    )
                    fig_tau_scatter1.update_layout(
                        plot_bgcolor="rgba(0, 0, 0, 0)",
                        paper_bgcolor="rgba(0, 0, 0, 0)",
                        xaxis_title="umap_1",
                        yaxis_title="umap_1",
                        title="All",
                    )
                    st.plotly_chart(
                        fig_tau_scatter1, use_container_width=True, height=1000
                    )
                with st_select_col2:
                    fig_tau_scatter2 = go.Figure(
                        data=build_scatter(
                            tau_adata,
                            color_umap_by=specificity_by,
                            pick_umap_color=[
                                display_tau_df.iloc[tau_select]["gene_name"].values[0]
                            ],
                        )
                    )
                    fig_tau_scatter2.update_layout(
                        plot_bgcolor="rgba(0, 0, 0, 0)",
                        paper_bgcolor="rgba(0, 0, 0, 0)",
                        xaxis_title="umap_1",
                        yaxis_title="umap_1",
                        title=f"{display_tau_df.iloc[tau_select]['gene_name'].values[0]} Expression",
                    )
                    st.plotly_chart(
                        fig_tau_scatter2, use_container_width=True, height=1000
                    )

        elif tau_or_label_specificity and "specificity" in tau_adata.uns.keys():
            specificity_df = tau_adata.uns["specificity"][specificity_by]
            st_col1, st_col2 = st.columns(2)
            with st_col1:
                fig_tau = go.Figure(
                    go.Histogram(
                        x=specificity_df.values.flatten()[::10],
                        nbinsx=50,
                    )
                )
                fig_tau.update_layout(
                    plot_bgcolor="rgba(0, 0, 0, 0)",
                    paper_bgcolor="rgba(0, 0, 0, 0)",
                    xaxis_title="specificity",
                )

                st.plotly_chart(fig_tau, use_container_width=True, height=1000)
                specificity_select = st.slider(
                    "Select Specificity:",
                    0.0,
                    1.0,
                    (0.25, 0.75),
                )
            with st_col2:
                display_selectivity_df = (
                    tau_adata.var.loc[
                        :,
                        ~tau_adata.var.columns.str.startswith("tau")
                        & (
                            ~tau_adata.var.columns.isin(
                                [
                                    "n_cells",
                                    "means",
                                    "dispersions",
                                    "dispersions_norm",
                                    "highly_variable",
                                    "feature_is_filtered",
                                    "feature_name",
                                    "feature_reference",
                                    "feature_biotype",
                                    "feature_length",
                                    "feature_type",
                                    "mean",
                                    "std",
                                ]
                            )
                        ),
                    ]
                    .assign(
                        specificity_score=specificity_df.max(axis=1),
                    )
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

                display_selectivity_df = display_selectivity_df.loc[
                    (
                        display_selectivity_df["specificity_score"]
                        >= specificity_select[0]
                    )
                    & (
                        display_selectivity_df["specificity_score"]
                        <= specificity_select[1]
                    )
                    & ~display_selectivity_df["specificity_score"].isna(),
                ].sort_values("specificity_score", ascending=False)
                st.write(
                    f"Genes: {display_selectivity_df.shape[0]}, uPs: {display_selectivity_df['vtx_id'].explode().shape[0]}"
                )
                specificity_select = st.dataframe(
                    display_selectivity_df,
                    use_container_width=True,
                    selection_mode="single-row",
                    on_select="rerun",
                )["selection"]["rows"]

            if specificity_select:
                st_select_col1, st_select_col2 = st.columns(2)
                with st_select_col1:
                    fig_specificity_scatter1 = go.Figure(
                        data=build_scatter(
                            tau_adata,
                            color_umap_by=specificity_by,
                            pick_umap_color=["All"],
                        )
                    )
                    fig_specificity_scatter1.update_layout(
                        plot_bgcolor="rgba(0, 0, 0, 0)",
                        paper_bgcolor="rgba(0, 0, 0, 0)",
                        xaxis_title="umap_1",
                        yaxis_title="umap_1",
                        title="All",
                    )
                    st.plotly_chart(
                        fig_specificity_scatter1, use_container_width=True, height=1000
                    )
                with st_select_col2:
                    fig_specificty_scatter2 = go.Figure(
                        data=build_scatter(
                            tau_adata,
                            color_umap_by=specificity_by,
                            pick_umap_color=[
                                display_selectivity_df.iloc[specificity_select][
                                    "gene_name"
                                ].values[0]
                            ],
                        )
                    )
                    fig_specificty_scatter2.update_layout(
                        plot_bgcolor="rgba(0, 0, 0, 0)",
                        paper_bgcolor="rgba(0, 0, 0, 0)",
                        xaxis_title="umap_1",
                        yaxis_title="umap_1",
                        title=f"{display_selectivity_df.iloc[specificity_select]['gene_name'].values[0]} Expression",
                    )
                    st.plotly_chart(
                        fig_specificty_scatter2, use_container_width=True, height=1000
                    )
        else:
            st.write("Under construction...:construction:")

        # with st.spinner("Thinking..."):
        #     specificity_metric, label = utils.label_expression_similarity(
        #         tau_adata.layers['ntpm'].copy(),
        #         tau_adata.obs.copy(),
        #         specificity_by,
        #     )

        # if "tau" in tau_adata.var.columns:
        #     fig_tau = go.Figure(
        #         go.Histogram(
        #             # x=tau_adata.var["tau"],
        #             x=specificity_metric.flatten(),
        #             nbinsx=25,
        #         )
        #     )
        #     fig_tau.update_layout(
        #         plot_bgcolor="rgba(0, 0, 0, 0)",
        #         paper_bgcolor="rgba(0, 0, 0, 0)",
        #         xaxis_title="tau",
        #     )

        #     st.plotly_chart(fig_tau, use_container_width=True, height=1000)
        #     tau_select = st.slider(
        #         "Select tau:",
        #         0.0,
        #         1.0,
        #         (0.25, 0.75),
        #     )
        #     tau_var = tau_adata.var.merge(
        #         velia_utils.pull_sorfensembl_mappings(),
        #         left_index=True,
        #         right_index=True,
        #         how="left",
        #     ).merge(
        #         velia_utils.pull_receptorensembl_mappings(),
        #         left_index=True,
        #         right_index=True,
        #         how="left",
        #     )
        #     display_tau_df = tau_var.loc[
        #         (tau_var["tau"] > tau_select[0])
        #         & (tau_var["tau"] < tau_select[1])
        #         & ~tau_var["tau"].isna(),
        #     ].sort_values("tau", ascending=False)
        #     st.write(
        #         f"Genes: {display_tau_df.shape[0]}, uPs: {display_tau_df['vtx_id'].explode().shape[0]}"
        #     )
        #     tau_select = st.dataframe(
        #         display_tau_df,
        #         use_container_width=True,
        #         selection_mode="single-row",
        #         on_select="rerun",
        #         column_config={
        #             "Select": st.column_config.CheckboxColumn(),
        #         },
        #     )["selection"]["rows"]
        #     st.write(display_tau_df.iloc[tau_select])

    # with st.expander("Marker Genes/Tissue Specificity"):
    #     choose_atlas_subset = st.radio(
    #         "Choose atlas subset:",
    #         ["Entire Atlas", "Selected Atlas Subset"],
    #         horizontal=True,
    #     )
    #     if choose_atlas_subset == "Entire Atlas":
    #         marker_adata = multiadata.full_atlas_adata
    #     else:
    #         marker_adata = adata

    #     marker_gene_by = st.radio(
    #         "Pick marker genes by:",
    #         ["Cell Type", "Tissue", "Tissue + Cell Type", "Louvain"],
    #         horizontal=True,
    #     )
    #     if marker_gene_by == "Tissue + Cell Type":
    #         marker_gene_by = "Tissue Cell Type"
    #     marker_gene_by = marker_gene_by.lower().replace(" ", "_")

    #     number_top_genes = int(
    #         st.number_input(
    #             "Number of top genes:",
    #             value=100,
    #             step=25,
    #             min_value=10,
    #             max_value=1000,
    #         )
    #     )

    #     pick_marker_cell = st.selectbox(
    #         f"Pick {marker_gene_by}:", marker_adata.obs[marker_gene_by].unique()
    #     )

    #     if (
    #         f"marker_selection_df_{choose_atlas_subset}_{marker_gene_by}_{resolution_option}"
    #         not in st.session_state
    #     ):
    #         with st.spinner("Thinking..."):
    #             marker_group_selection_df = (
    #                 pd.concat(
    #                     [
    #                         marker_adata.obs.copy(),
    #                         pd.DataFrame(
    #                             marker_adata.layers["ntpm"].copy(),
    #                             index=marker_adata.obs.index.copy(),
    #                             columns=marker_adata.var["gene_name"].copy(),
    #                         ),
    #                     ],
    #                     axis=1,
    #                     ignore_index=False,
    #                 )
    #                 .groupby(marker_gene_by)[marker_adata.var["gene_name"]]
    #                 .agg("mean")
    #                 .fillna(0.0)
    #             )
    #             marker_group_selection_df = (
    #                 marker_group_selection_df - marker_group_selection_df.mean(axis=0)
    #             ) / marker_group_selection_df.std(axis=0)
    #         st.session_state[
    #             f"marker_selection_df_{choose_atlas_subset}_{marker_gene_by}_{resolution_option}"
    #         ] = marker_group_selection_df
    #     else:
    #         marker_group_selection_df = st.session_state[
    #             f"marker_selection_df_{choose_atlas_subset}_{marker_gene_by}_{resolution_option}"
    #         ]

    #     top_genes = marker_group_selection_df.loc[pick_marker_cell, :].sort_values(
    #         ascending=False
    #     )[:number_top_genes]
    #     heatmap_df = marker_group_selection_df.loc[:, top_genes.index]
    #     fig_ts_heatmap = go.Figure(
    #         data=go.Heatmap(
    #             z=heatmap_df.values,
    #             x=top_genes.index,
    #             y=heatmap_df.index,
    #             colorbar_title="standardized nTPM",
    #         ),
    #     )
    #     fig_ts_heatmap.update_layout(
    #         xaxis_tickmode="linear",
    #         yaxis_tickmode="linear",
    #     )
    #     st.plotly_chart(fig_ts_heatmap, use_container_width=True, height=2000)
