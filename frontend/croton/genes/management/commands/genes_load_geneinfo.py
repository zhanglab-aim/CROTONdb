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

from genes.models import Gene, CrossRefDB, CrossRef, Alias
from organisms.models import Organism


ENTREZ_URL = 'ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/'
ORG = 'Homo_sapiens'
SUFFIX = 'gene_info.gz'


class Command(BaseCommand):
    def add_arguments(self,parser):
        parser.add_argument('--geneinfo_file',help='Input filname')
        parser.add_argument('--remote', action="store_true",help = 'Download gene_info file remotely from NCBI Entrez')
        parser.add_argument('--tax_id', help='Input filname',type=int)
        parser.add_argument('--gi_tax_id',help = "To work with cerivisiae's tax id change, use, this, otherwise we will just use tax_id")
        parser.add_argument('--symbol_col',  help = "The column containing the symbol id.", default = 2,type=int)
        parser.add_argument('--systematic_col',type=int, help = "The column containing the systematic id.  If this is '-' or blank, the symbol will be used", default = 3)
        parser.add_argument('--alias_col',type=int,help = "The column containing gene aliases.  If this is '-' or blank, the symbol will be used", default = 4)

    # option_list = BaseCommand.option_list + (
    #     make_option('--geneinfo_file', action = 'store', dest = 'geneinfo'),
    #     make_option('--remote', action = 'store_true', dest = 'remote', help = 'Download gene_info file remotely from NCBI Entrez', default=False),
    #     make_option('--tax_id', action = 'store', dest = 'tax_id'),
    #     make_option('--gi_tax_id', action = 'store', dest = 'gi_tax_id', help = "To work with cerivisiae's tax id change, use, this, otherwise we will just use tax_id"),
    #     make_option('--symbol_col', action = 'store', dest = 'symbol_col', help = "The column containing the symbol id.", default = 2),
    #     make_option('--systematic_col', action = 'store', dest = 'systematic_col', help = "The column containing the systematic id.  If this is '-' or blank, the symbol will be used", default = 3),
    #     make_option('--alias_col', action = 'store', dest = 'alias_col', help = "The column containing gene aliases.  If this is '-' or blank, the symbol will be used", default = 4),
    # )
    #
    help = 'Add standards from stds_file with the associations from assoc_file.'
    def handle(self, *args, **options):
        print(options)
        tax_id = options['tax_id']
        geneinfo = options['geneinfo_file']
        remote = options['remote']
        gi_tax_id = tax_id
        symb_col = options['symbol_col']
        syst_col = options['systematic_col']
        alias_col = options['alias_col']
        if options['gi_tax_id']:
            gi_tax_id = ['gi_tax_id']

        if geneinfo:
            geneinfo = open(geneinfo)
        elif remote:
            print(ENTREZ_URL + ORG + '.' + SUFFIX)
            gene_zip_fh = urllib.request.urlopen(ENTREZ_URL + ORG + '.' + SUFFIX, timeout=5)
            print(gene_zip_fh)
            geneinfo = gzip.GzipFile(fileobj=io.BytesIO(gene_zip_fh.read()))
            print(geneinfo)
            gene_zip_fh.close()
            print("closed")
        print("about to load")
        print(tax_id)
        print(geneinfo)
        if tax_id and geneinfo:
            gene_info = {}
            org_matches = 0
            for line in geneinfo:
                toks = line.strip().split('\t')
                if toks[0] == str(gi_tax_id):
                    org_matches += 1
                    (entrezid, symbol, systematic_name, aliases, crossrefs, description, status) = (toks[1], toks[symb_col], toks[syst_col], toks[alias_col], toks[5], toks[8], toks[9])
                    if status == 'protein-coding' or status == 'rRNA':
                        #ignoring any pseudogenes without a symbol
                        if status == 'pseudo' and symbol.startswith('LOC'):
                            continue
                        if status == 'hypothetical protein':
                            continue
                        if (not systematic_name) or (systematic_name == '-'):
                            systematic_name = symbol
                        if (not entrezid) or (entrezid == '-'):
                            continue
                        else:
                            try:
                                gene = gene_info[entrezid]
                            except KeyError:
                                gene = gene_info[entrezid] = {}
                        gene['systematic_name'] = systematic_name
                        if symbol and symbol != '-':
                            gene['standard_name'] = symbol
                        if aliases and (aliases != '-'):
                            try:
                                alias_set = gene['aliases']
                            except KeyError:
                                alias_set = gene['aliases'] = set()
                            alias_set.update([str(x) for x in aliases.split('|')])
                        if crossrefs and (crossrefs != '-'):
                            try:
                                xref_set = gene['xrefs']
                            except KeyError:
                                xref_set = gene['xrefs'] = set()
                            xref_set.update([str(x) for x in crossrefs.split('|')])
                        if description and (description != '-'):
                            try:
                                other_desc = gene['description']
                                gene['description'] = other_desc + ' and ' + description
                            except KeyError:
                                gene['description'] = description

            non_unique_crossrefs = []
            org = Organism.objects.get(taxid = tax_id)
            gene_set = set(Gene.objects.filter(organism = org).values_list('entrez', flat = True))
            xrdbcache = {}
            bum_xrdbs = set()
            logger.info("About to save %d genes" % len(gene_info.keys()))
            for entrez in gene_info.keys():
                gdict = gene_info[entrez]
                try:
                    gene = Gene.objects.get(entrez=entrez, organism=org)
                    if (gene.standard_name != gdict['standard_name']) or (gene.description != gdict['description']) or (gene.systematic_name != gdict['systematic_name']):
                        gene.systematic_name = gdict['systematic_name']
                        gene.standard_name = gdict['standard_name']
                        gene.description = gdict['description']
                        gene.save()
                    if gene.entrez in gene_set:
                        gene_set.remove(gene.entrez)
                    else:
                        logger.error('Gene with systematic_name %s not in gene_set', gene.systematic_name)
                except Gene.DoesNotExist:
                    logger.debug('Gene %s did not exist and was created for %s', gdict['standard_name'], tax_id)
                    gene = Gene(entrez=entrez, systematic_name=gdict['systematic_name'], standard_name=gdict['standard_name'], description=gdict['description'], organism=org)
                    try:
                        with transaction.atomic():
                            gene.save()
                    except DatabaseError:
                        logger.error('%s had some trouble with the gene_info record.', systematic_name)

                try:
                    new_xrefs = gdict['xrefs']
                except KeyError:
                    new_xrefs = set()
                    logger.info('%s had no xrefs', entrez)

                cur_xrefs = CrossRef.objects.filter(gene = gene).select_related('crossrefdb')
                cur_xref_proc = set()
                for item in cur_xrefs:
                    if (item.crossrefdb.name + ':' + item.xrid) not in new_xrefs:
                        logger.info('The cross reference %s from %s was removed for %s', item.xrid, item.crossrefdb.name, gdict['standard_name'])
                        item.delete()
                    else:
                        cur_xref_proc.add(item.crossrefdb.name + ':' + item.xrid)
                if new_xrefs is not None:
                    for item in new_xrefs.difference(cur_xref_proc):
                        (xrname, xrid) = item.split(':',1)
                        if xrname in bum_xrdbs:
                            continue
                        try:
                            xrdb = xrdbcache[xrname]
                        except KeyError:
                            try:
                                xrdb = CrossRefDB.objects.get(name = xrname)
                            except CrossRefDB.DoesNotExist:
                                bum_xrdbs.add(xrname)
                                logger.error('CrossRefDB %s did not exist', xrname, extra = {'xrname': xrname})
                                continue
                            xrdbcache[xrname] = xrdb

                        try:
                            xref = CrossRef.objects.get(crossrefdb = xrdb, gene = gene, xrid = xrid)
                        except CrossRef.DoesNotExist:
                            try:
                                with transaction.atomic():
                                    xref = CrossRef(crossrefdb = xrdb, gene = gene, xrid = xrid)
                                    logger.debug('The cross reference %s in %s was added for %s', item, xrname, gdict['standard_name'])
                                    xref.save()
                            except IntegrityError:
                                non_unique_crossrefs.append((xrdb, xrid))
                try:
                    alias = Alias.objects.get(gene = gene, organism = org)
                    if 'aliases' not in gdict:
                        alias.delete()
                    else:
                        new_aliases_txt = ','.join(gdict['aliases'])
                        if alias.name != new_aliases_txt:
                            alias.name = new_aliases_txt
                            alias.save()
                except Alias.DoesNotExist:
                    if 'aliases' in gdict:
                        try:
                            logger.debug('Aliases for %s did not exist and was created for %s', gdict['standard_name'], tax_id)
                            alias = Alias(organism = org, gene = gene, name = ','.join(list(gdict['aliases'])))
                            with transaction.atomic():
                                alias.save()
                        except DatabaseError:
                            logger.error('%s had some trouble with the gene_info record.', systematic_name)

            for gene in gene_set:
                Gene.objects.filter(entrez=gene, organism = org).delete()
                logger.warning('The gene %s no longer appears in the gene_info file and was removed.', gene)
            for (xrdb, xrid) in non_unique_crossrefs:
                CrossRef.objects.filter(crossrefdb = xrdb, xrid = xrid).delete()
                logger.info('The cross reference %s in %s was not unique and was removed.', xrid, xrdb)
            geneinfo.close()
        else:
            logger.error('Couldn\'t load geneinfo %s for org %s.', options.get('geneinfo'), tax_id , exc_info = sys.exc_info(), extra = {'options': options})




