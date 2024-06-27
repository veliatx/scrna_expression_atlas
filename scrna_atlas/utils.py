from pathlib import Path
import shutil
from typing import *
from urllib.request import urlretrieve
import logging
import warnings

import numpy as np
import anndata as ad
import scanpy as sc
import pandas as pd
import pyranges as pr
from s3fs import S3FileSystem
from scipy import sparse
from tqdm import tqdm
import alphashape
import streamlit as st

from shapely.ops import unary_union
from shapely import Polygon, MultiPolygon, GeometryCollection

from scrna_atlas import settings

logger = logging.getLogger(__name__)


def move_cellxgene_adatas(
    cellxgene_loc: Path,
    dest_loc: Path,
    title_substring: str,
) -> None:
    """ """

    for al in cellxgene_loc.glob("*.h5ad"):
        adata = ad.read_h5ad(al, backed="r")
        if title_substring in adata.uns["title"]:
            shutil.copy(
                al,
                dest_loc / cellxgene_loc.name,
            )


def download_cellxgene_adatas(
    dest_loc: Path = settings.DATA_LOC,
    url: str = settings.CELLXGENE_ADATA_URL,
    check_exists: bool = True,
) -> None:
    """ """
    pass


def download_hpa_scrna_data(
    dest_loc: Path = settings.HPA_DATA_LOC,
    check_exists: bool = True,
) -> None:
    """ """
    hpa_scrna_url = Path(settings.HPA_SCRNA_URL)
    if check_exists and (dest_loc / hpa_scrna_url.name).exists():
        logger.info(f"Already downloaded {hpa_scrna_url.name}.")
    else:
        try:
            urlretrieve(hpa_scrna_url, dest_loc / hpa_scrna_url.name)
        except Exception as e:
            logger.exception(e)
            raise Exception(f"Failed to download {hpa_scrna_url}.")

    hpa_scrna_cluster_url = Path(settings.HPA_SCRNA_CLUSTER_URL)
    if check_exists and (dest_loc / hpa_scrna_cluster_url.name).exists():
        logger.info(f"Already downloaded {hpa_scrna_cluster_url.name}.")
    else:
        try:
            urlretrieve(hpa_scrna_cluster_url, dest_loc / hpa_scrna_cluster_url.name)
        except Exception as e:
            logger.exception(e)
            raise Exception(f"Failed to download {hpa_scrna_cluster_url}.")

    hpa_scrna_cluster_expression_url = Path(settings.HPA_SCRNA_CLUSTER_EXPRESSION_URL)
    if check_exists and (dest_loc / hpa_scrna_cluster_expression_url.name).exists():
        logger.info(f"Already downloaded {hpa_scrna_cluster_expression_url.name}.")
    else:
        try:
            urlretrieve(
                hpa_scrna_cluster_expression_url,
                dest_loc / hpa_scrna_cluster_expression_url.name,
            )
        except Exception as e:
            logger.exception(e)
            raise Exception(f"Failed to download {hpa_scrna_cluster_expression_url}.")


def pull_genes_gtf(
    gtf_fh: Union[Path, str],
) -> pd.DataFrame:
    """ """
    df = pr.read_gtf(gtf_fh).df

    df = df.loc[:, ["gene_id", "gene_name", "gene_type"]]
    df["gene_id_noversion"] = df["gene_id"].str.split(".").str[0]
    df.set_index(["gene_id_noversion"], inplace=True)
    df = df.loc[~df.index.duplicated(keep="first")]

    return df


def tmm_normalize(
    data: np.ndarray,
    trim_lfc=0.3,
    trim_mag=0.05,
) -> Tuple[np.ndarray, np.ndarray]:
    """TMM-normalization, reworked from conorm package.

    Args:
        data (np.ndarray): Expression matrix, rows as samples and columns as genes.
        trim_lfc (float): Cutoff for fold change.
        trim_mag (float): Cutoff for magnitude.
    Returns:
        (Tuple[np.ndarray,np.ndarray]): Normalized expression values, and normalization factors.
    """
    x = data.copy()

    lib_size = np.nansum(x, axis=1)
    mask = x == 0
    x[:, mask.all(axis=0)] = np.nan
    p75 = np.nanpercentile(x, 75, axis=1)
    ref = np.argmin(np.abs(p75 - p75.mean()))
    mask[:, mask[ref]] = True
    x[mask] = np.nan

    norm_x = x / lib_size[:, np.newaxis]
    log = np.log2(norm_x)
    m_g = log - log[ref]
    a_g = (log + log[ref]) / 2

    perc_m_g = np.nanquantile(
        m_g,
        [trim_lfc, 1 - trim_lfc],
        axis=1,
        method="nearest",
    )[..., np.newaxis]

    perc_a_g = np.nanquantile(
        a_g,
        [trim_mag, 1 - trim_mag],
        axis=1,
        method="nearest",
    )[..., np.newaxis]

    mask = mask | (m_g < perc_m_g[0]) | (m_g > perc_m_g[1])
    mask = mask | (a_g < perc_a_g[0]) | (a_g > perc_a_g[1])
    w_gk = (1 - norm_x) / x
    w_gk = 1 / (w_gk + w_gk[ref])

    w_gk[mask] = 0
    m_g[mask] = 0
    w_gk = w_gk / np.nansum(w_gk, axis=1)[:, np.newaxis]
    tmms = np.nansum(w_gk * m_g, axis=1)
    tmms = tmms - tmms.mean()
    tmms = np.exp2(tmms)

    return data / tmms[..., np.newaxis], tmms


def collapse_counts_simpleaf(
    adata: ad.AnnData,
    method: str = "snRNA",
    remove_versions: bool = True,
    sparsify: bool = True,
    make_unique_varnames: bool = True,
) -> ad.AnnData:
    """ """
    splice_locs = ~adata.var.index.str.endswith(("U", "A"))
    unsplice_locs = adata.var.index.str.endswith("U")
    ambiguous_locs = adata.var.index.str.endswith("A")

    adata.var["gene_id"] = adata.var.index.str.split("-").str[0]

    adata_splice = adata[:, splice_locs].copy()
    adata_unsplice = adata[:, unsplice_locs].copy()
    adata_ambiguous = adata[:, ambiguous_locs].copy()

    adata_splice.var.index = adata_splice.var["gene_id"]
    adata_unsplice.var.index = adata_unsplice.var["gene_id"]
    adata_ambiguous.var.index = adata_ambiguous.var["gene_id"]

    var = pd.DataFrame(
        adata.var["gene_id"].unique(),
        columns=["gene_id"],
        index=adata.var["gene_id"].unique(),
    )

    if remove_versions:
        var.index = var.index.str.split(".").str[0]

    adata_out = ad.AnnData(
        X=sparse.csr_matrix(
            (adata.shape[0], adata.var["gene_id"].nunique()), dtype=np.float32
        ),
        obs=adata.obs.copy(),
        var=var.copy(),
    )
    if not np.array_equal(
        adata_splice.var["gene_id"], adata_unsplice.var["gene_id"]
    ) and np.array_equal(adata_unsplice.var["gene_id"], adata_ambiguous.var["gene_id"]):
        raise Exception(
            "Relative ordering of splice, unpsliced, and ambiguous transcripts needs to be the same."
        )

    if method == "snRNA":
        adata_out.X = adata_splice.X + adata_unsplice.X + adata_ambiguous.X
    elif method == "scRNA":
        adata_out.X = adata_splice.X + adata_ambiguous.X
    else:
        raise Exception("Not a valid method.")

    if isinstance(adata_out.X, np.ndarray) and sparsify:
        adata_out.X = sparse.csr_matrix(adata_out.X)
        adata_out.X.eliminate_zeros()

    adata_out = adata_out[:, adata_out.X.sum(axis=0) != 0.0].copy()

    if make_unique_varnames:
        adata_out.var_names_make_unique()

    return adata_out


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

    if isinstance(polygon, MultiPolygon) or isinstance(polygon, GeometryCollection):
        for g in polygon.geoms:
            x.extend([p[0] for p in g.exterior.coords] + [None])
            y.extend([p[1] for p in g.exterior.coords] + [None])
    else:
        x.extend([p[0] for p in polygon.exterior.coords] + [None])
        y.extend([p[1] for p in polygon.exterior.coords] + [None])

    return x, y


def convert_scatter_alphashape(
    pts: np.ndarray,
    alpha: float = 10.0,
) -> None:
    """ """
    a = alphashape.alphashape(
        pts,
        alpha,
    )
    if isinstance(a, MultiPolygon) or isinstance(a, GeometryCollection):
        if len(list(a.geoms)) == 0:
            return np.array([])
        a = sorted(a.geoms, key=lambda x: x.area)[-1]

    return np.vstack(
        [
            a.exterior.xy[0],
            a.exterior.xy[1],
        ]
    )


def pseudobulk_adata(
    adata: ad.AnnData,
    bulking_key: str,
    celltype_key: str,
    pseudoreplicate_key: Union[str, None] = None,
    pseudoreplicate_n: int = 1,
    min_cells: int = 10,
    min_counts: int = 1000,
    aggregate_method: str = "sum",
    umap_key: str = "X_scvi_umap",
    seed: int = 42,
    polygon_alpha: float = 5.0,
    calc_tpm: bool = True,
    normalize_length: bool = True,
) -> ad.AnnData:
    """ """
    if aggregate_method not in ("sum", "mean", "median"):
        raise Exception("aggregate_method must be one of: sum, mean, median.")

    if pseudoreplicate_key and pseudoreplicate_n > 1:
        if pseudoreplicate_key not in adata.obs.columns:
            raise KeyError(
                f"Need to specify a valid column in obs as pseudoreplicate_key: {pseudoreplicate_key}."
            )

        pseudoreplicate_df = adata.obs.copy()
        pseudoreplicate_df = pseudoreplicate_df.sample(frac=1, random_state=seed)

        pseudoreplicate_df["pseudoreplicate"] = pseudoreplicate_df.groupby(
            pseudoreplicate_key
        ).cumcount()
        pseudoreplicate_df["pseudoreplicate"] = pseudoreplicate_df.apply(
            lambda x: f'{x[pseudoreplicate_key]}__{x["pseudoreplicate"] % pseudoreplicate_n}',
            axis=1,
        )

        adata.obs = adata.obs.merge(
            pseudoreplicate_df.loc[:, ["pseudoreplicate"]],
            left_index=True,
            right_index=True,
        )

        adata.obs["_sample_key"] = adata.obs.apply(
            lambda x: f'{x[celltype_key]}__{x[bulking_key]}{"__"+x["pseudoreplicate"]}',
            axis=1,
        )
    else:
        adata.obs["_sample_key"] = adata.obs.apply(
            lambda x: f"{x[celltype_key]}__{x[bulking_key]}",
            axis=1,
        )

    adata.obs["_cell_index"] = np.arange(len(adata.obs))

    pobs = adata.obs.groupby("_sample_key").agg("first").reset_index()
    pobs.set_index("_sample_key", inplace=True)

    for g, d in tqdm(
        adata.obs.groupby("_sample_key"), total=adata.obs["_sample_key"].nunique()
    ):
        cell_indices = d["_cell_index"].values
        cells = adata.X[cell_indices, :]

        pobs.loc[g, "bulk_n_cells"] = cells.shape[0]
        pobs.loc[g, "bulk_total_cells"] = cells.sum()

    pobs = pobs.loc[
        (pobs["bulk_n_cells"] >= min_cells) & (pobs["bulk_total_cells"] >= min_counts)
    ]

    pX = np.zeros((len(pobs), adata.X.shape[1]))
    puns = {
        "title": adata.uns["title"] if "title" in adata.uns.keys() else None,
        "citation": adata.uns["citation"] if "citation" in adata.uns.keys() else None,
        "assay": (
            ",".join(adata.obs["assay"].unique())
            if "assay" in adata.obs.columns
            else None
        ),
        "umap_polygons": {},
    }
    players = {}
    players["percent_expression"] = np.zeros_like(pX)
    players["count_mean"] = np.zeros_like(pX)
    players["count_std"] = np.zeros_like(pX)

    for i, (n, r) in tqdm(enumerate(pobs.iterrows()), total=pobs.shape[0]):
        indices = adata.obs.loc[adata.obs["_sample_key"] == n, "_cell_index"]
        players["percent_expression"][i] = np.divide(
            (adata.X[indices, :] > 0).sum(axis=0), indices.size
        )
        players["count_mean"][i] = adata.X[indices, :].toarray().mean(axis=0)
        players["count_std"][i] = adata.X[indices, :].toarray().std(axis=0)

        if aggregate_method == "sum":
            pX[i] = adata.X[indices, :].sum(axis=0)
        elif aggregate_method == "mean":
            pX[i] = adata.X[indices, :].mean(axis=0)
        elif aggregate_method == "median":
            pX[i] = np.median(adata.X[indices, :], axis=0)
        puns["umap_polygons"][n] = convert_scatter_alphashape(
            adata.obsm[umap_key][n == adata.obs["_sample_key"], :],
            polygon_alpha,
        )

    if calc_tpm:
        if normalize_length:
            samples_sums = np.nansum(
                (pX / np.array(adata.var["feature_length"])[np.newaxis, ...]), axis=1
            )
            players["tpm"] = (
                1e6
                * (pX / np.array(adata.var["feature_length"])[np.newaxis, ...])
                / samples_sums[..., np.newaxis]
            )
            players["ntpm"], pobs["tpm_tmm_norm_factor"] = tmm_normalize(
                players["tpm"].copy()
            )
        else:
            sample_sums = np.nansum(pX, axis=1)
            players["tpm"] = 1e6 * pX / sample_sums[..., np.newaxis]
            players["ntpm"], pobs["tpm_tmm_norm_factor"] = tmm_normalize(
                players["tpm"].copy()
            )

    pobs.drop("_cell_index", inplace=True, axis=1)

    padata = ad.AnnData(
        pX.copy(),
        obs=pobs.copy(),
        var=adata.var.copy(),
        layers=players.copy(),
        uns=puns.copy(),
    )

    return padata


def csv_to_csr_matrix(
    fh: Path,
) -> Tuple[sparse.csr_matrix, List[str], List[str]]:
    """Helper function for reading HPA single-cell, otherwise these blow up in
    memory because they're stored as dense arrays."""

    genes = []
    data = []
    indices = []
    indptr = [0]

    with open(fh, "r") as f_in:
        cell_ids = f_in.readline().strip().split("\t")[1:]
        for i, r in enumerate(f_in):
            line = r.strip().split("\t")
            genes.append(line[0])
            row = [
                (ri, float(c))
                for ri, c in enumerate(line[1:])
                if not c.startswith(("0", "0.0"))
            ]
            data.extend([r[1] for r in row])
            indices.extend([r[0] for r in row])
            indptr.append(indptr[-1] + len(row))

    m = sparse.csr_matrix(
        (data, indices, indptr),
        dtype=int,
        shape=(
            len(genes),
            len(cell_ids),
        ),
    )

    return m.transpose(), genes, cell_ids


def calc_adata_tau(
    adata: ad.AnnData,
    layer: str = "ntpm",
    tau_key: Union[str, None] = None,
    return_tau: bool = False,
    log_transform: bool = True,
) -> None:
    """ """
    if log_transform:
        l = np.log2(adata.layers[layer] + 1.0)
    else:
        l = adata.layers[layer].copy()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        xi = l / l.max(axis=0)
        tau = (1 - xi).sum(axis=0) / (xi.shape[0] - 1)
    if not tau_key:
        tau_key = "tau"
    if not return_tau:
        adata.var[tau_key] = tau
    else:
        return tau


def label_expression_similarity(
    adata_layer: np.ndarray,
    adata_obs: pd.DataFrame,
    adata_var: pd.DataFrame,
    specificity_column: Union[str, None] = None,
) -> Tuple[np.ndarray, np.ndarray]:
    """ """
    dummy_specifity_df = pd.get_dummies(
        adata_obs[specificity_column] if specificity_column else adata_obs.index
    ).astype(float)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        exp_norm = np.log2(adata_layer + 1) / np.linalg.norm(
            np.log2(adata_layer + 1), axis=0
        )
        exp_norm[~np.isfinite(exp_norm)] = 0.0
        label_norm = dummy_specifity_df.values / np.linalg.norm(
            dummy_specifity_df.values, axis=0
        )
        s = np.dot(exp_norm.T, label_norm)

    specificity_df = pd.DataFrame(
        s, columns=dummy_specifity_df.columns, index=adata_var.index.copy()
    )
    return specificity_df


def calculate_hpastyle_enrichments(
    adata: ad.AnnData,
    num_top: int = 5,
    fc_thresh: float = 4,
    obs_column: Union[None, str] = None,
) -> None:
    """ """
    for i, (ix, r) in enumerate(adata.var.iterrows()):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            exp = adata.layers["ntpm"][:, i] / adata.layers["ntpm"][:, i].max()
            sorted_exp = sorted(
                adata.layers["ntpm"][:, i] / adata.layers["ntpm"][:, i].max(),
                reverse=True,
            )
            tissue_enriched = sorted_exp[1] <= (1 / fc_thresh)
            for n in list(range(2, num_top + 1))[::-1]:
                group_enriched = (
                    (sorted_exp[n:] / sorted_exp[n - 1]) <= (1 / fc_thresh)
                ).all()
                group_enrich_tissue = (
                    adata.obs.index[exp.argsort()[::-1][:n]]
                    if not obs_column
                    else adata.obs[obs_column][exp.argsort()[::-1][:n]]
                )
                if group_enriched:
                    break
            tissue_enhanced_mask = (
                adata.layers["ntpm"][:, i].mean() / adata.layers["ntpm"][:, i]
            ) <= (1 / fc_thresh)
            tissue_enhanced = tissue_enhanced_mask.any()
        tissue_enhanced_tissue = (
            adata.obs.index[tissue_enhanced_mask]
            if not obs_column
            else adata.obs[obs_column][tissue_enhanced_mask]
        )
        tissue_enriched_tissue = group_enrich_tissue[0]
        adata.var.loc[
            ix,
            [
                "enriched",
                "group_enriched",
                "enhanced",
                "not_specific",
                "enriched_tissue",
                "group_enriched_tissue",
                "tissue_enhanced_tissue",
            ],
        ] = [
            tissue_enriched,
            (~tissue_enriched & group_enriched),
            (~tissue_enriched & ~group_enriched & tissue_enhanced),
            (~tissue_enriched & ~group_enriched & ~tissue_enhanced),
            tissue_enriched_tissue if tissue_enriched else np.NaN,
            (
                ",".join(group_enrich_tissue)
                if (~tissue_enriched & group_enriched)
                else np.NaN
            ),
            (
                ",".join(tissue_enhanced_tissue)
                if (~tissue_enriched & ~group_enriched & tissue_enhanced)
                else np.NaN
            ),
        ]


def process_cellxgene_atlas_adatas(
    adata: ad.AnnData,
    adata_name: str,
    out_loc: Path,
    use_rep: str,
    umap_key: str,
    feature_lengths_df: Union[None, pd.DataFrame] = None,
    normalize_length: bool = False,
) -> None:
    """ """
    sc.pp.neighbors(
        adata,
        n_neighbors=15,
        n_pcs=50,
        use_rep=use_rep,
    )

    sc.tl.louvain(
        adata,
        resolution=1.0,
    )

    sc.tl.umap(adata)

    adata.obs["tissue_cell_type"] = adata.obs.apply(
        lambda x: f"{x.tissue}__{x.cell_type}".replace(" ", "_"),
        axis=1,
    )

    if "feature_length" not in adata.var.columns:
        print(adata.var.shape)
        adata.var = adata.var.merge(
            feature_lengths_df,
            left_index=True,
            right_index=True,
            how="left",
        )
        print(adata.var.shape)

    padata = pseudobulk_adata(
        adata,
        "louvain",
        "tissue_cell_type",
        min_cells=5,
        min_counts=200,
        umap_key=umap_key,
        calc_tpm=True,
        normalize_length=normalize_length,
    )

    padata.var.loc[:, "gene_name"] = padata.var.loc[:, "feature_name"]

    padata.write_h5ad(out_loc / f"{adata_name}_pseuodobulk.h5ad")
