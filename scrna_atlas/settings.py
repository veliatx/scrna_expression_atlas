from pathlib import Path

import numpy as np

DATA_LOC = Path("/home/ec2-user/scrna/data")
GENOME_LOC = Path("/efs/expression_atlas/scrna/genomes")
GTF_HUMAN_FH = GENOME_LOC / "hg38" / "gencode.v38.primary_assembly.annotation.gtf"
GTF_MOUSE_FH = GENOME_LOC / "M27" / "gencode.vM27.primary_assembly.annotation.gtf"
EXPERIMENT_LOC = DATA_LOC / "experiments"
EXPERIMENT_EXPLORER_INDEX = EXPERIMENT_LOC / "experiments.tsv"
CELLXGENE_ADATA_URL = ""
HPA_SCRNA_URL = "https://www.proteinatlas.org/download/rna_single_cell_read_count.zip"
HPA_SCRNA_CLUSTER_URL = (
    "https://www.proteinatlas.org/download/rna_single_cell_cluster_description.tsv.zip"
)
HPA_SCRNA_CLUSTER_EXPRESSION_URL = (
    "https://www.proteinatlas.org/download/rna_single_cell_type_tissue.tsv.zip"
)
HPA_DATA_LOC = DATA_LOC / "hpa"

DASHBOARD_SORFS_URL = "s3://velia-data-dev/VDC_005_veliadb_backup/dashboard_data/20240610_v_1_11/cache/sorf_df.parq"
RECEPTORS_URL = "s3://velia-analyses-dev/VAP_20240306_receptor_coverage/data/receptors_uniprot_ccle_expression.csv"

DE_LFC_THRESH = np.log2(1.5)
DE_PADJ_THRESH = 0.05
