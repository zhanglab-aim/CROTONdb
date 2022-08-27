import logging
import sys

logger = logging.getLogger(__name__)

from django.core.management.base import BaseCommand
from genes.models import CrossRefDB, CrossRef, Gene


class Command(BaseCommand):
    def add_arguments(self,parser):
        parser.add_argument('--uniprot_file',help = 'filtered uniprot file (i.e. with: zgrep "GeneID" idmapping.dat.gz > uniprot_entrez.txt'),

    help = 'Add uniprot identifiers to database.'
    def handle(self, *args, **options):
        uniprot_file = options['uniprot_file']
        if uniprot_file:
            uniprot_file = open(uniprot_file)
        else:
            logger.error("Please provide uniprot file --uniprot_file")
            sys.exit(0)

        if uniprot_file:
            entrez_set = set(Gene.objects.values_list('entrez', flat = True))
            logger.info("Total entrez ids: %s", len(entrez_set))
            uniprot_entrez = {}
            for line in uniprot_file:
                (uniprot_id, itype, entrez_id) = line.strip().split()
                if int(entrez_id) in entrez_set:
                    uniprot_entrez[uniprot_id] = entrez_id
            uniprot = CrossRefDB.objects.get(name = 'UniProtKB')
            new_xrs = []
            logger.info("Total uniprot ids: %s", len(uniprot_entrez.keys()))
            for uniprot_id in uniprot_entrez.keys():
                try:
                    uniprot_xr = CrossRef.objects.get(crossrefdb = uniprot, xrid = uniprot_id)
                    uniprot_xr.gene_id = uniprot_entrez[uniprot_id]
                    uniprot_xr.save()
                except CrossRef.DoesNotExist:
                    new_xrs.append(CrossRef(crossrefdb=uniprot, xrid=uniprot_id, gene_id=uniprot_entrez[uniprot_id]))
            CrossRef.objects.bulk_create(new_xrs)
            uniprot_file.close()




