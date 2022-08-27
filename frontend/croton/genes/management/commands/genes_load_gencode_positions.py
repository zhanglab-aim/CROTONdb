import logging
import sys
logger = logging.getLogger(__name__)
logger.addHandler(logging.StreamHandler())

import urllib.request
import gzip
import io

from django.core.management.base import BaseCommand
from django.db import IntegrityError, transaction
from django.db.utils import DatabaseError

from genes.models import Gene, CrossRefDB, CrossRef, GeneInterval


## Input file is encode.<v>.basic.annnotation.gtf
## E.g. http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_35/gencode.v35.basic.annotation.gtf.gz


def parse_ensemble_id(input_str):
    """
    input_str - string like this: gene_id "ENSG00000187608.10"
    output  - ENSG00000187608"""
    ensemble_id_str = input_str.strip().split(' ')[1]
    ensemble_id = (ensemble_id_str.strip('"')).split('.')[0]
    return ensemble_id

def parse_gene_name(input_str):
    """
    input_str - string like this: gene_name "ISG15"
    output  - ISG15
    """
    gene_name = (input_str.strip().split(' ')[1]).strip('"')
    return gene_name

class Command(BaseCommand):
    def add_arguments(self,parser):
        parser.add_argument('--input',help='Input filname',required=True)

    help = 'Add gencode gene positions'
    def handle(self, *args, **options):
        geneinfo = options['input']

        gi_count = GeneInterval.objects.count()
        if gi_count > 0:
            print("WARNING! GeneInterval table is not empty. There are {} rows in the table.".format(gi_count))
            val = input(" Is it OK to delete all data? [y/n] ")
            if val != 'y':
                print("Exiting")
                return
            GeneInterval.objects.all().delete()
            print(" done")


        print("\nQuerying all ensemble ids")
        #query all EnsembleIds
        ensemble_gene_dict = {}
        # for xr in CrossRef.objects.filter(crossrefdb__name='Ensembl'):
        #     ensemble_gene_dict[xr.xrid] = xr.gene

        xrid_list = ['ENSG00000284770','ENSG00000285053']
        #for xr in CrossRef.objects.filter(xrid__in=xrid_list):
        for xr in CrossRef.objects.filter(crossrefdb__name='Ensembl'):        
            ensemble_gene_dict[xr.xrid] = xr.gene


        print("Querying all entrez ids")
        gene_dict = {}
        #for gene in Gene.objects.filter(entrez__in=[6905,219899,116804918]):
        for gene in Gene.objects.all():    
            gene_dict[gene.entrez] = gene
        print("\tdone")

        geneinfo = open(geneinfo)
        #to store the information by chromosomes
        chrom_gene_info = {}        
        for line in geneinfo:
            #skip comments
            if line.startswith('##'):
                continue
            #['chr1', 'HAVANA', 'gene', '11869', '14409', '.', '+', '.', 'gene_id "ENSG00000223972.5"; gene_type "transcribed_unprocessed_pseudogene"; gene_name "DDX11L1"; level 2; hgnc_id "HGNC:37102"; havana_gene "OTTHUMG00000000961.2";']
            toks = line.strip().split('\t')                
            chrom = toks[0]
            annot_type = toks[2]
            if not annot_type == 'gene':
                continue
            start = int(toks[3])
            end = int(toks[4])
            annot_str = toks[8]
            annots = annot_str.split(';')            
            if not annots[1].strip() == 'gene_type "protein_coding"':
                continue
            gene_id_annot = annots[0].strip()
            gene_name_annot = annots[2].strip()    
            ensemble_id = parse_ensemble_id(gene_id_annot)     
            gene_name = parse_gene_name(gene_name_annot)  
            #print("\t{} {}\n".format(gene_id_annot,gene_name_annot))
            #print("\t{} {}\n".format(ensemble_id,gene_name))
            if not chrom in chrom_gene_info:
                chrom_gene_info[chrom] = {}
            chrom_gene_info[chrom][ensemble_id] = {}
            chrom_gene_info[chrom][ensemble_id]["chrom"]  = chrom
            chrom_gene_info[chrom][ensemble_id]["name"] = gene_name
            chrom_gene_info[chrom][ensemble_id]["start"] = start
            chrom_gene_info[chrom][ensemble_id]["end"] = end
        
        geneinfo.close()

        print("All lines read, there are {} chromosomes found".format(len(chrom_gene_info.keys())))

        

        #genes that have more than 1 ensemble ids 
        duplicate_entries = {}
        #Saving has to be done by chromosome because if there are genes that have multiple start/end positions
        #they need to be combined
        for chrom in chrom_gene_info.keys():       
            gene_interval_dict = {}
            gene_info = chrom_gene_info[chrom]            
            print("\nThere are {} gene positions found for chromosome >{}<".format(len(gene_info.keys()),chrom))
            #print(gene_info)
            for ensemble in gene_info.keys():
                gdict = gene_info[ensemble]    
                gene = None
                if ensemble in ensemble_gene_dict:
                    gene = ensemble_gene_dict[ensemble]
                    #print("Gene found {} for {}".format(gene,ensemble))                    
                    # gene_intervals_list.append(GeneInterval(gene=gene,chromosome=gene_info[ensemble]["chrom"],start=gene_info[ensemble]["start"],
                    #                                         end=gene_info[ensemble]["end"]))

                    
                else:
                    gene_name = gene_info[ensemble]["name"]
                    print("\tMissing gene-ensemble id for {} ({})".format(ensemble,gene_name))
                    try:
                        gene = Gene.objects.get(standard_name=gene_name)
                    except Gene.DoesNotExist:
                        continue
                    print("\t\tfound: {}({})".format(gene.standard_name,gene.entrez))


                if not gene.entrez in gene_interval_dict:
                    gene_interval_dict[gene.entrez] = {}
                    gene_interval_dict[gene.entrez]['start_pos'] = []
                    gene_interval_dict[gene.entrez]['end_pos'] = []
                gi_dict = gene_interval_dict[gene.entrez]
                gi_dict['start_pos'].append(gene_info[ensemble]["start"])
                gi_dict['end_pos'].append(gene_info[ensemble]["end"])



            gene_intervals_list = []
            #time to save the gene_intervals
            for entrez in gene_interval_dict:
                gi_dict = gene_interval_dict[entrez]

                #if there are more than 1 start/end positions, take the minimum start and maximum start position
                start_position = min(gi_dict['start_pos'])
                end_position = max(gi_dict['end_pos'])

                gene = gene_dict[entrez]

                if(len(gi_dict['start_pos']) > 1):
                    duplicate_entries['{}-{}'.format(entrez,gene.standard_name)] = gi_dict

                gene_intervals_list.append(GeneInterval(gene=gene,chromosome=chrom,start=start_position,end = end_position))
            
            print('Saving {} gene intervals for {}'.format(len(gene_intervals_list),chrom))
            
            #print(gene_intervals_list)
            GeneInterval.objects.bulk_create(gene_intervals_list)


        if len(duplicate_entries.keys()) > 0:
            print('Duplicate entries')
            for entry in duplicate_entries.keys():
                print('{}: {}'.format(entry,duplicate_entries[entry]))

