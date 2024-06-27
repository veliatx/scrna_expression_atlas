from pathlib import Path

import pandas as pd
import anndata as ad


class AtlasAdataWrapper:

    def __init__(
        self,
        adatas_locs: Path,
        atlas_name: str,
        full_atlas: str = None,
        glob_pattern: str = "*.h5ad",
    ):
        """
        Args:
            adatas_loc (Path): Path to folder with all of the adatas.
            atlas_name (str): Display name of atlas.
            full_atlas (str): Name of the full atlas, used for picking out the full atlas adata.
            glob_pattern (str): Glob pattern used for globbing adatas out of folder.
        """
        self._atlas_name = atlas_name
        self._adatas_locs = adatas_locs
        self._full_atlas_adata = full_atlas
        self.adatas = {}
        self._glob_pattern = glob_pattern
        self.load_adatas()

    @property
    def list_adatas(self):
        return [
            k
            for k in self.adatas.keys()
            if "-ALL" in k.replace(" ", "").upper()
            or "ALL-" in k.replace(" ", "").upper()
        ] + [
            k
            for k in self.adatas.keys()
            if not "-ALL" in k.replace(" ", "").upper()
            and not "ALL-" in k.replace(" ", "").upper()
        ]

    @property
    def full_atlas_adata(self):
        if self._full_atlas_adata:
            return self.adatas[self._full_atlas_adata]
        else:
            return self.adatas[self.list_adatas[0]]

    def load_adatas(self):
        """TODO: Figure out if there is a way to open these with a context manager so we
        don't have a bunch of handles to the adata files. Ideally would like to read all in
        off disk to avoid lots of memory overhead."""
        for al in self._adatas_locs.glob(self._glob_pattern):
            adata = ad.read_h5ad(al, backed="r")
            self.adatas[adata.uns["title"]] = adata
