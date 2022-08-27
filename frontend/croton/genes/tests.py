from django.test.client import Client
from django.test import TestCase
from django.core import management

from genes.models import Gene, CrossRefDB
from organisms.models import Organism

class GeneTest(TestCase):
    def setUp(self):
        self.human = Organism.objects.create(taxid=9606, name="Homo sapiens")
        self.gene = Gene.objects.create(entrez=1, systematic_name="A1BG", 
                standard_name="A1BG", description="alpha-1-B glycoprotein", 
                organism=self.human)

    def tearDown(self):
        Organism.objects.all().delete()
        Gene.objects.all().delete()


    def test_gene(self):
        """
        Tests that save creates the organism's slug.
        """
        gene = Gene.objects.get(entrez=1)
        self.assertEqual(gene.standard_name, "A1BG")    

    def test_management(self):
        """
        Tests the management commands that load genes.
        """
        management.call_command('genes_add_xrdb', name="Ensembl", URL="http://www.ensembl.org/Gene/Summary?g=_REPL_")
        management.call_command('genes_add_xrdb', name="Vega", URL="http://vega.sanger.ac.uk/Homo_sapiens/Gene/Summary?g=_REPL_")
        management.call_command('genes_add_xrdb', name="Entrez", URL="http://www.ncbi.nlm.nih.gov/gene/_REPL_")
        management.call_command('genes_add_xrdb', name="HGNC", URL="http://www.genenames.org/data/hgnc_data.php?hgnc_id=_REPL_")
        management.call_command('genes_add_xrdb', name="HPRD", URL="http://www.hprd.org/protein/_REPL_")
        management.call_command('genes_add_xrdb', name="UniProtKB", URL="http://www.uniprot.org/uniprot/_REPL_")
        management.call_command('genes_add_xrdb', name="MIM", URL="http://www.omim.org/entry/_REPL_")
        #management.call_command('genes_update_from_uniprot', uniprot_file="test/uniprot_entrez.test")

        CrossRefDB.objects.get(name='Ensembl')

        management.call_command('genes_load_geneinfo', geneinfo_file="genes/test/human_gene_info", tax_id='9606', symbol_col=2, systematic_col=2)

        self.assertTrue(Gene.objects.all().count() == 8)


    def test_duplicate_geneinfo(self):
    #     """
    #     Tests the management commands that load genes.
    #     """
        management.call_command('genes_add_xrdb', name="Ensembl", URL="http://www.ensembl.org/Gene/Summary?g=_REPL_")
        management.call_command('genes_add_xrdb', name="Vega", URL="http://vega.sanger.ac.uk/Homo_sapiens/Gene/Summary?g=_REPL_")
        management.call_command('genes_add_xrdb', name="Entrez", URL="http://www.ncbi.nlm.nih.gov/gene/_REPL_")
        management.call_command('genes_add_xrdb', name="HGNC", URL="http://www.genenames.org/data/hgnc_data.php?hgnc_id=_REPL_")
        management.call_command('genes_add_xrdb', name="HPRD", URL="http://www.hprd.org/protein/_REPL_")
        management.call_command('genes_add_xrdb', name="UniProtKB", URL="http://www.uniprot.org/uniprot/_REPL_")
        management.call_command('genes_add_xrdb', name="MIM", URL="http://www.omim.org/entry/_REPL_")
    
        management.call_command('genes_load_geneinfo', geneinfo_file="genes/test/human_gene_info_bad", tax_id='9606', symbol_col=2, systematic_col=2)

        self.assertTrue(Gene.objects.all().count() == 3)


