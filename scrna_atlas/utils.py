from pathlib import Path
import shutil

import anndata as ad
from s3fs import S3FileSystem

from scrna_atlas import settings

def move_cellxgene_adatas(
    cellxgene_loc: Path,
    dest_loc: Path, 
    title_substring:str,
) -> None:
    """ """

    for al in cellxgene_loc.glob('*.h5ad'):
        adata = ad.read_h5ad(al, backed='r')
        if title_substring in adata.uns['title']:
            shutil.copy(
                al,
                dest_loc / cellxgene_loc.name,
            )

def download_cellxgene_adatas(
    dest_loc: Path=settings.DATA_LOC,
    url: str=settings.CELLXGENE_ADATA_URL,
) -> None:
    """ """
    pass

def download_hpa_scrna_data(
    
) -> None:
    """ """
