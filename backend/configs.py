import os


LOCAL_DIR = os.path.dirname(os.path.realpath(__file__))
DATA_DIR = f"{LOCAL_DIR}/../frontend/data"

GENOME_FA_PATH = f"{DATA_DIR}/genomes/GRCh38.primary_assembly.genome.fa"
GTF_DB_PATH = f"{DATA_DIR}/genomes/gencode.v35.annotation.gff3.gz.gtf_sqldb"

# GENOME_FA_PATH = f"{DATA_DIR}/genomes/Small_Test_Genome.all.fa"
# GTF_DB_PATH = f"{DATA_DIR}/genomes/Small_Test_Genome.CROTONTest.gff3.gtf_sqldb"


CDS_PAM_DIR = f"{DATA_DIR}/bed/CDSpams-byChrom/"

GNOMAD_AF_ENTRIES = ['AF', 'AF_male', 'AF_female', 'AF_afr', 'AF_afr_female', 'AF_afr_male', 'AF_amr', 'AF_amr_female', 'AF_amr_male',
    'AF_asj', 'AF_asj_female', 'AF_asj_male', 'AF_eas', 'AF_eas_female', 'AF_eas_jpn', 'AF_eas_kor', 'AF_eas_male', 'AF_eas_oea',
    'AF_fin', 'AF_fin_female', 'AF_fin_male', 'AF_nfe', 'AF_nfe_bgr', 'AF_nfe_est', 'AF_nfe_female',
    'AF_nfe_male', 'AF_nfe_nwe', 'AF_nfe_onf',  'AF_nfe_seu', 'AF_nfe_swe', 'AF_oth', 'AF_oth_female', 'AF_oth_male',
    'AF_popmax', 'AF_raw', 'AF_sas', 'AF_sas_female', 'AF_sas_male']


GENE_BED_DIR = f"{DATA_DIR}/bed/byGene/"

