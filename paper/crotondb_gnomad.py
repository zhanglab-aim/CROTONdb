import pandas as pd
from tqdm import tqdm
from gnomad_db.database import gnomAD_DB # https://github.com/KalinNonchev/gnomAD_DB/tree/master

data_folder = '/home/zhangz2/zhanglab/fzzhang/CROTONdb/old_datas/052121-variant-2'
database_location = "test_dir"
db = gnomAD_DB(database_location, gnomad_version="v3")

def query_maf(res):
    return var_df2

chroms = [str(_) for _ in range(1,23)] + ['X']
#chroms = ['4']
res = []
for chrom in chroms:
    data_fp = f'{data_folder}/tabix/{chrom}/CROTON_varpred_{chrom}.gz'
    print(f"processing {chrom}..")
    try:
        with pd.read_table(data_fp, chunksize=10**5) as reader:
            for df in tqdm(reader):
                res.append(df.query('diff_frameshift>0.2 or diff_frameshift<-0.2'))
    except KeyboardInterrupt:
        break
    var_df = pd.concat(res)
    var_df.reset_index(drop=True, inplace=True)
    var_df['chrom'] = var_df['chrom'].str.replace('chr', '')
    var_df['pos'] = var_df['pos'].astype(int)
    var_df.to_csv(f'crotondb_out/{chrom}.nofreq.csv')

    freqs = db.get_info_from_df(var_df[['chrom', 'pos', 'ref', 'alt']], "AF_eas, AF_nfe, AF_fin, AF_afr, AF_asj")
    freqs['pop_max'] = freqs.max(axis=1)
    freqs['pop_min'] = freqs.min(axis=1)

    var_df2 = pd.concat([var_df, freqs], axis=1)
    var_df2.to_csv(f'crotondb_out/{chrom}.freq.csv')



