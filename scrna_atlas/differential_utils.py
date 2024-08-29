from pathlib import Path
from itertools import combinations
import logging
from typing import *

import numpy as np
import anndata as ad
import pandas as pd
from s3fs import S3FileSystem
from scipy import sparse
from tqdm import tqdm
from joblib import Parallel, delayed, parallel_backend

from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from pydeseq2.utils import build_design_matrix, nb_nll

from rpy2.robjects.packages import importr
import rpy2.robjects as ro

from scrna_atlas import settings

logger = logging.getLogger(__name__)

limma = importr("limma")
limma.contrastfit = ro.r("contrasts.fit")
edger = importr("edgeR")
deseq2 = importr("DESeq2")


def run_pseudobulk_de(
    padata: ad.AnnData,
    factors: List[str],
    ref_factor: Dict[str, str],
    groupby: str = "cell_type",
    groups: Union[List[str], None] = None,
    frac_expressing: float = 1 / 3,
    min_counts: int = 5,
    cell_factor_threshold: int = 2,
    method: str = "pydeseq2",
) -> None:
    """ """
    if not isinstance(groups, list):
        groups = []

    for ct in padata.obs[groupby].unique():

        sub_adata = padata[padata.obs[groupby] == ct].copy()
        sub_adata = sub_adata[
            :,
            (sub_adata.X.sum(axis=0) > min_counts)
            & (
                (sub_adata.X > min_counts).sum(axis=0)
                > int(np.round(sub_adata.shape[0] * frac_expressing))
            ),
        ].copy()
        cell_factor_counts = sub_adata.obs[[groupby] + factors].value_counts()
        if not (cell_factor_counts > cell_factor_threshold).all():
            logger.info(f"Skipping {ct} because not enough samples/cell type.")
            continue

        if not all([r in sub_adata.obs[c].tolist() for c, r in ref_factor.items()]):
            logger.info(f"Skipping {ct} because reference factor not present in {ct}.")
            continue

        if method.lower() == "pydeseq2":
            dds = pydeseq2_fit_model(
                sub_adata,
                factors,
                ref_factor,
                groups,
            )
            pydeseq2_test_contrasts(
                padata,
                dds,
                factors,
                ref_factor,
                de_prefix=ct,
            )
        elif method.lower() == "deseq2":
            model, _ = r_build_model(
                sub_adata,
                factors,
                ref_factor,
                groups=groups,
            )
            r_dds = r_deseq2_fit_model(
                sub_adata,
                model,
                factors,
                groups,
            )
            r_deseq2_test_contrasts(
                padata,
                model,
                r_dds,
                de_prefix=ct,
            )
            for k, df in padata.varm.items():
                if not k.startswith("de-"):
                    continue
                df["log10_padj"] = np.log10(df["padj"])

        elif method.lower() == "edger":
            model, _ = r_build_model(
                sub_adata,
                factors,
                ref_factor,
                groups=groups,
            )
            r_dge = r_edger_fit_model(
                sub_adata,
                model,
            )
            r_edger_test_contrasts(
                padata,
                model,
                r_dge,
                de_prefix=ct,
            )
            for k, df in padata.varm.items():
                if not k.startswith("de-"):
                    continue
                df.rename(
                    {
                        "logFC": "log2FoldChange",
                        "logCPM": "baseMean",
                        "F": "stat",
                        "PValue": "pvalue",
                        "FDR": "padj",
                    },
                    axis=1,
                    inplace=True,
                )
                df["log10_padj"] = np.log10(df["padj"])
        else:
            raise Exception("Method must be one of {pydeseq2, deseq2, edger}.")


def pydeseq2_fit_model(
    padata: ad.AnnData,
    factors: List[str],
    ref_factor: Dict[str, str],
    groups: Union[List[str], None] = None,
) -> DeseqDataSet:
    """ """
    ref_level = [list(ref_factor.keys())[0], ref_factor[list(ref_factor.keys())[0]]]
    dds = DeseqDataSet(
        adata=padata,
        design_factors=factors + (groups if groups else []),
        ref_level=ref_level,
        refit_cooks=True,
    )
    dds.deseq2()

    return dds


def pydeseq2_test_contrasts(
    padata: ad.AnnData,
    dds: DeseqDataSet,
    factors: List[str],
    ref_factor: Dict[str, str],
    de_prefix: str = "",
) -> None:
    """ """

    for f in factors:
        for cl in dds.obs[f].unique():
            if cl == ref_factor[f]:
                continue
            sub_ds = DeseqStats(dds, contrast=[f, cl, ref_factor[f]])
            sub_ds.summary()
            sub_ds.lfc_shrink()

            padata.varm[
                f"de-{de_prefix}__{cl}_vs_{ref_factor[f]}".replace(" ", "")
            ] = padata.var.loc[:, ["gene_name"]].merge(
                sub_ds.results_df,
                left_index=True,
                right_index=True,
                how="left",
            )


def relevel_design(
    dds: DeseqDataSet,
    ref_level: Tuple[str, str],
) -> None:
    """Relevels pydeseq2 DeseqDataSet to level in ref_level. Rearranges coefficients to accomodate
    new reference level and rebuilds design matrix.

    Args:
        dds (DeseqDataset): pydeseq2 object to modify.
        ref_level (Tuple[str,str]): Tuple of two elements ('condition','new_reference_level').
    """
    if ref_level[0] not in dds.obs.columns:
        raise ValueError("%s condition not in design." % ref_level[0])
    if ref_level[1] not in dds.obs[ref_level[0]].values:
        raise ValueError("%s condition level not in %s." % (ref_level[1], ref_level[0]))
    if not dds.ref_level:
        raise AttributeError(
            '%s define reference level "ref_level" for original design.'
        )

    if any(
        True if c.startswith(ref_level[0]) and c.endswith(ref_level[1]) else False
        for c in dds.obsm["design_matrix"].columns
    ):
        logger.info("%s already reference level for %s" % (ref_level[0], ref_level[1]))
        return

    design_matrix = build_design_matrix(
        metadata=dds.obs.copy(),
        design_factors=dds.design_factors,
        ref_level=ref_level,
    )

    if "LFC" not in dds.varm.keys():
        dds.deseq2()

    coef_df = dds.varm["LFC"].copy()

    refo = [c.split("_")[-1] for c in dds.varm["LFC"] if c.startswith(ref_level[0])][0]

    columns_to_relevel = [
        c for c in design_matrix.columns if c.startswith(ref_level[0])
    ]

    coef_df["intercept"] = (
        dds.varm["LFC"]["intercept"]
        + dds.varm["LFC"]["%s_%s_vs_%s" % (ref_level[0], ref_level[1], refo)]
    )

    for c in columns_to_relevel:
        ref, con = c.split(ref_level[0] + "_")[-1].split("_vs_")

        if "%s_%s_vs_%s" % (ref_level[0], con, ref) in dds.varm["LFC"].columns:
            coef_df[c] = (
                -1.0 * dds.varm["LFC"]["%s_%s_vs_%s" % (ref_level[0], con, ref)]
            )
        else:
            coef_df[c] = (
                dds.varm["LFC"]["%s_%s_vs_%s" % (ref_level[0], ref, refo)]
                - dds.varm["LFC"]["%s_%s_vs_%s" % (ref_level[0], con, refo)]
            )

    columns_drop = [c for c in dds.varm["LFC"] if c.startswith(ref_level[0])]

    coef_df.drop(columns_drop, axis=1, inplace=True)
    coef_df = coef_df[design_matrix.columns]

    dds.varm["LFC"] = coef_df.copy()
    dds.obsm["design_matrix"] = design_matrix.copy()
    dds.ref_level = ref_level

    logger.info("dds releveled to %s-%s" % (ref_level[0], ref_level[1]))


def r_build_model(
    padata: ad.AnnData,
    factors: List[str],
    ref_factor: Dict[str, str],
    groups: Union[List[str], None] = None,
    factorial_design: bool = False,
) -> Tuple[ro.vectors.FloatMatrix, Dict[str, List[str]]]:
    """
    Args:
        padata (ad.AnnData): The pseudobulk anndata.
        factors (List[str]): List of columns in padata.obs to assemble design matrix.
        ref_factor (Dict[str,str]):
        groups (Union[List[str], None]):
    Returns:
    """

    for f, r in ref_factor.items():
        if f not in padata.obs.columns.tolist():
            raise KeyError(f"Reference factor {f} not in padata.obs.columns.")

        if r not in padata.obs[f].tolist():
            raise ValueError(f"Reference factor {r} not a factor in padata.obs[{f}]")

    for f in factors:
        if f not in padata.obs.columns.tolist():
            raise KeyError(f"Factor {f} not in padata.obs.columns.")

    if groups:
        for g in groups:
            if g not in padata.obs.columns.tolist():
                raise KeyError(f"Group {f} not in padata.obs.columns.")

    # Create factor vectors to pass to model.matrix.

    factor_vectors = {}
    for f in factors:
        if f in ref_factor.keys():
            # Set reference factor as the first in the factor vector.
            factor_levels = [ref_factor[f]] + [
                l
                for i, l in enumerate(padata.obs[f])
                if l not in padata.obs[f][:i].tolist() and l != ref_factor[f]
            ]
        else:
            factor_levels = [
                l
                for i, l in enumerate(padata.obs[f])
                if l not in padata.obs[f][:i].tolist()
            ]
        factor_vector = ro.FactorVector(
            padata.obs[f].tolist(), levels=ro.StrVector(factor_levels)
        )
        factor_vectors[f] = factor_vector

    if groups:
        # Create group vectors to pass to model.matrix.
        group_vectors = {}
        for g in groups:
            group_levels = [
                l for i, l in enumerate(padata.obs[g]) if l not in padata.obs[g][:i]
            ]

            group_vector = ro.FactorVector(
                padata.obs[g].tolist(), levels=ro.StrVector(group_levels)
            )
            group_vectors[g] = group_vector

    # Create model formula to generate design matrix.
    if not factorial_design:
        model_formula = "~" + "-".join(factors)
    else:
        model_formula = (
            "~"
            + "+".join(factors)
            + "+"
            + "+".join([f"{c1}:{c2}" for c1, c2 in combinations(factors, 2)])
        )

    if groups:
        model_formula += "+" + "+".join(groups)

    fmla = ro.Formula(model_formula)
    for f, v in factor_vectors.items():
        fmla.environment[f] = v

    if groups:
        for g, v in group_vectors.items():
            fmla.environment[g] = v

    model = ro.r("model.matrix")(fmla)

    # Rename the columns of the design matrix.
    colnames = ["Intercept"]
    for f, v in factor_vectors.items():
        levels = v.levels
        colnames.extend([f"{l}_vs_{levels[0]}" for l in levels[1:]])

    if factorial_design:
        for c1, c2 in combinations(factors, 2):
            levels_1 = factor_vectors[c1].levels
            levels_2 = factor_vectors[c2].levels

            # In R, the second condition in interaction term is iterated over first,
            # then the first condition in the interaction with the first condition displayed first.

            for l2 in levels_2[1:]:
                for l1 in levels_1[1:]:
                    colnames.append(f"{l1}_vs_{levels_1[0]}:{l2}_vs_{levels_2[0]}")

    if groups:
        for g, v in group_vectors.items():
            levels = v.levels
            colnames.extend([f"group{l}__{levels[0]}" for l in levels[1:]])

    # Replace columns with correct column names.
    model.colnames = ro.StrVector(colnames)

    # Convert factor_vectors back to python object.
    for k in factor_vectors.keys():
        factor_vectors[k] = list(factor_vectors[k].levels)

    return model, factor_vectors


def r_build_expression_matrix(
    padata: ad.AnnData,
    layer: Union[str, None] = None,
    use_raw: bool = False,
    obsm: Union[str, None] = None,
) -> ro.vectors.FloatMatrix:
    """ """
    if layer:
        if layer not in padata.layers.keys():
            raise KeyError(f"Layer {layer} not found in padata.layers.")
        if isinstance(padata.layers[layer], np.ndarray):
            X = padata.layers[layer].astype(float).copy()
        else:
            try:
                X = padata.layers[layer].to_array().astype(float)
            except Exception as e:
                raise Exception("Cannot convert layer to dense.") from e
    if obsm:
        if obsm not in padata.obsm.keys():
            raise KeyError(f"Key {obsm} not found in padata.obsm")
        try:
            X = padata.obsm[obsm].astype(float).copy()
        except Exception as e:
            raise Exception("Cannot convert obsm key.") from e
    elif use_raw:
        if isinstance(padata.raw.X, np.ndarray):
            X = padata.raw.X.astype(float).copy()
        else:
            try:
                X = padata.raw.X.to_array().astype(float)
            except Exception as e:
                raise Exception("Cannot conver X to dense.") from e
    else:
        if isinstance(padata.X, np.ndarray):
            X = padata.X.astype(float).copy()
        else:
            try:
                X = padata.X.to_array().astype(float)
            except Exception as e:
                raise Exception("Cannot conver X to dense.") from e

    r_m = ro.r.matrix(
        ro.FloatVector(X.T.flatten().tolist()),
        ncol=padata.obs.shape[0],
    )

    r_m.colnames = ro.StrVector(padata.obs.index.values)
    r_m.rownames = ro.StrVector(padata.var.index.values)

    return r_m


def r_limma_fit_model(
    padata: ad.AnnData,
    model: ro.vectors.FloatMatrix,
    layer: Union[str, None] = None,
) -> ro.vectors.ListVector:
    """ """
    r_m = r_build_expression_matrix(
        padata,
        layer=layer,
    )

    fit = limma.lm_fit(r_m, model)

    return fit


def r_limma_test_contrasts(
    padata: ad.AnnData,
    model: ro.vectors.FloatMatrix,
    r_fit: ro.vectors.ListVector,
    de_prefix: str = "",
    varm_prefix: str = "de",
) -> None:
    """ """
    for i, c in enumerate(model.colnames):
        if c == "Intercept":
            continue
        r_cfit = limma.contrastsfit(r_fit, i)
        r_cfit = limma.eBayes(r_cfit)
        r_tt = limma.topTable(r_cfit, n="Inf")

        res_df = pd.DataFrame(
            np.array(r_tt).T,
            index=r_tt.rownames,
            columns=r_tt.colnames,
        )

        padata.varm[f"{varm_prefix}-{de_prefix}__{c}"] = padata.var.loc[
            :, ["gene_name"]
        ].merge(
            res_df,
            left_index=True,
            right_index=True,
            how="left",
        )


def r_edger_fit_model(
    padata: ad.AnnData,
    model: ro.vectors.FloatMatrix,
    layer: Union[str, None] = None,
) -> ro.vectors.ListVector:
    """ """
    r_m = r_build_expression_matrix(
        padata,
        layer=layer,
    )

    r_dge = edger.DGEList(r_m)
    r_dge = edger.calcNormFactors(r_dge)
    r_dge = edger.estimateDisp(r_dge, design=model)

    fit = edger.glmQLFit(r_dge, model)

    return fit


def r_edger_test_contrasts(
    padata: ad.AnnData,
    model: ro.vectors.FloatMatrix,
    r_fit: ro.vectors.ListVector,
    de_prefix: str = "",
) -> None:
    """ """
    for i, c in enumerate(model.colnames):
        if c == "Intercept":
            continue
        r_qlf = edger.glmQLFTest(r_fit, coef=i + 1)
        r_tt = edger.topTags(r_qlf, n="Inf")

        res_df = pd.DataFrame(
            np.array(r_tt[0]).T,
            index=r_tt[0].rownames,
            columns=r_tt[0].colnames,
        )

        padata.varm[f"de-{de_prefix}__{c}"] = padata.var.loc[:, ["gene_name"]].merge(
            res_df,
            left_index=True,
            right_index=True,
            how="left",
        )


def r_deseq2_fit_model(
    padata: ad.AnnData,
    model: ro.vectors.FloatMatrix,
    factors: List[str],
    groups: List[str],
    layer: Union[str, None] = None,
) -> ro.methods.RS4:
    """ """
    r_m = r_build_expression_matrix(
        padata,
        layer=layer,
    )
    # Not sure how necessary creating a full colData is if we're
    # passing in a design matrix.
    r_meta = ro.r.matrix(
        ro.StrVector(padata.obs.loc[:, factors + groups].values.flatten().tolist()),
        ncol=len(factors + groups),
    )
    r_meta.rownames = ro.StrVector(padata.obs.index.astype(str))
    r_meta.colnames = ro.StrVector(factors + groups)
    r_dds = deseq2.DESeqDataSetFromMatrix(
        countData=r_m,
        colData=r_meta,
        design=model,
    )
    r_dds = deseq2.DESeq(r_dds)

    return r_dds


def r_deseq2_test_contrasts(
    padata: ad.AnnData,
    model: ro.vectors.FloatMatrix,
    r_dds: ro.methods.RS4,
    de_prefix: str = "",
) -> None:
    """ """

    for i, c in enumerate(model.colnames):
        if c == "Intercept":
            continue
        r_res = deseq2.lfcShrink(r_dds, coef=i + 1, type="apeglm")
        # r_dds.do_slot_assign("condition", ro.r("relevel")(r_dds.do_slot("condition"), model.colnames[-1]))
        # r_res = deseq2.results(r_dds, contrast=[model.colnames[-1], model.colnames[-2]], )# type='normal')

        # Could probably push this through rpy2 pandas converter, but
        # have had difficulties getting these conversions to go well in the past.

        res_df = pd.DataFrame(
            np.array(r_res.do_slot("listData")).T,
            index=np.array(r_res.do_slot("rownames")),
            columns=np.array(r_res.do_slot("listData").do_slot("names")),
        )

        padata.varm[f"de-{de_prefix}__{c}"] = padata.var.merge(
            res_df,
            left_index=True,
            right_index=True,
            how="left",
        )
