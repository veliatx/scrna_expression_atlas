from pathlib import Path

import pandas as pd
import anndata as ad


class AtlasAdataWrapper:

    def __init__(self, adatas_locs: Path, atlas_name: str):
        self._atlas_name = atlas_name
        self._adatas_locs = adatas_locs
        self.adatas = {}
        self.load_adatas()

    @property
    def list_adatas(self):
        return [
            k for k in self.adatas.keys() if '- ALL' in k.upper() or 'ALL - ' in k.upper()
        ] + [
            k for k in self.adatas.keys() if not '- ALL' in k.upper() and not 'ALL - ' in k.upper()
        ]

    def load_adatas(self):
        for al in self._adatas_locs.glob('*.h5ad'):
            adata = ad.read_h5ad(al, backed='r')
            self.adatas[adata.uns['title']] = adata


