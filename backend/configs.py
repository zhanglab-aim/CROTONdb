import os

LOCAL_DIR = os.path.dirname(os.path.realpath(__file__))
ROOT_DIR = f"{LOCAL_DIR}/.."
DATA_DIR = f"{LOCAL_DIR}/../frontend/data"

# Fasta file location (*.fa)
GENOME_FA_PATH = f"{DATA_DIR}/genomes/GRCh38.primary_assembly.genome.fa"

# GTF DB file location (*.gff3.gz.gtf_sqldb)
GTF_DB_PATH = f"/common/zhangz2lab/shared/genomes/gff/gencode.v35.annotation.gff3.gz.gtf_sqldb"

# Bed file location (*.bed and *.intersect.bed)
CDS_PAM_DIR = f"{DATA_DIR}/bed/CDSpams-byChrom/"

# (separate by gene name) bed file location
GENE_BED_DIR = f"{DATA_DIR}/bed/byGene/"

# model location (CROTON.h5)
MODEL_PATH = os.path.join(ROOT_DIR, "CROTON.h5")

# tsv dir, store the output file with model prediction
TSV_DIR = os.path.join(ROOT_DIR, 'datavl', 'variant', 'tsv')

# tabix dir, store the tabix file location
TABIX_DIR = os.path.join(ROOT_DIR, 'datavl', 'variant', 'tabix')

GNOMAD_AF_ENTRIES = ['AF', 'AF_male', 'AF_female', 'AF_afr', 'AF_afr_female', 'AF_afr_male', 'AF_amr', 'AF_amr_female', 'AF_amr_male',
    'AF_asj', 'AF_asj_female', 'AF_asj_male', 'AF_eas', 'AF_eas_female', 'AF_eas_jpn', 'AF_eas_kor', 'AF_eas_male', 'AF_eas_oea',
    'AF_fin', 'AF_fin_female', 'AF_fin_male', 'AF_nfe', 'AF_nfe_bgr', 'AF_nfe_est', 'AF_nfe_female',
    'AF_nfe_male', 'AF_nfe_nwe', 'AF_nfe_onf',  'AF_nfe_seu', 'AF_nfe_swe', 'AF_oth', 'AF_oth_female', 'AF_oth_male',
    'AF_popmax', 'AF_raw', 'AF_sas', 'AF_sas_female', 'AF_sas_male']
