# Read Me


## Layout
The folder `data-builder` contains jupyter notebooks for compiling test-scale datsets.

The folder `slurm` contains slurm scripts for submitting batch jobs to HPC.

These scripts will use some source data stored in `frontend/data` folder. See `configs.py` for more detailed data filepaths.
You can soft-link your existing files to these folders.


## Workflow

1. get all PAM sites in coding regions (CDS)

2. intersect PAM to a target VCF file

3. Run ML predictors on reference and alternate sequences
    - CROTON
    - FORECasT
    - InDelphi
    - SPROUT

4. Store in a database

5. Serve in Django



## Downloading External Data Files

| Filename                                           | URL                                                                                                 | Notes                                |
|----------------------------------------------------|-----------------------------------------------------------------------------------------------------|--------------------------------------|
| GRCh38.primary_assembly.genome.fa                  | https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_35/GRCh38.p13.genome.fa.gz        | hg38 Fasta genome sequences          |
| gencode.v35.annotation.gff3.gz.gtf                 | https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_35/gencode.v35.annotation.gff3.gz | Gencode V35 annotation in GFF format |
| gnomad.exomes.r2.1.1.sites.liftover_grch38.vcf.bgz | https://gnomad.broadinstitute.org/downloads#v2-liftover-variants                                    | Gnomad v2 Exome variants; rename from .bgz to .gz             |