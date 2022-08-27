import re

from django.db import models
from django.urls import reverse
from organisms.models import Organism

"""
Genes will be added from an online database like 'Entrez', using the Django management
commands (see genes/management/commands/genes_load_geneinfo.py).
"""

class Gene(models.Model):
    # NCBI gene entrez id is used as the primary key
    entrez = models.PositiveIntegerField(primary_key=True)

    # Used for organisms like yeast
    systematic_name = models.CharField(max_length=32, db_index=True)

    # Gene symbol will be standard_name or systematic_name
    # if there is no standard_name
    standard_name = models.CharField(max_length=32, null=True)

    description = models.TextField()
    organism = models.ForeignKey(Organism, on_delete=models.CASCADE)

    # Gene ids from other bioinformatic databases
    xrefs = []

    # Gene aliases provided in the entrez gene_info file
    aliases = []

    class Meta:
        ordering = ['organism', 'systematic_name']

    def get_absolute_url(self):
        return reverse('genes:detail', args=[self.entrez])

    def __str__(self):
        if self.standard_name:
            return self.standard_name
        else:
            return self.systematic_name

    def get_xrefs(self):
        if not self.xrefs:
            self.xrefs = self.crossref_set.select_related('crossrefdb').all()
        return self.xrefs

    
    def get_gene_intervals(self):
        return self.geneinterval_set.all()

    def get_aliases(self):
        if not self.aliases:
            try:
                self.aliases = Alias.objects.get(gene=self)
                self.aliases = self.aliases.name.split(',')
            except Alias.DoesNotExist:
                pass

        return self.aliases

    def as_dict(self, more={}):
        xrefs = [x.as_dict() for x in self.get_xrefs()]
        gdict = {'standard_name': self.standard_name, 'systematic_name': self.systematic_name,
                 'entrez': self.entrez, 'description': self.description,
                 'xrefs': xrefs, 'aliases': self.get_aliases()}
        gdict.update(more)
        return gdict

    def wall_of_name(self):
        '''
        Appends identifiers for the different databases (such as Entrez id's) and returns them.
        Uses the CrossRef class below.
        '''
        names = [str(self.entrez), self.standard_name, self.systematic_name]
        names.extend([xref.xrid for xref in self.crossref_set.all()])
        for i in range(len(names)):
            names[i] = re.sub(r'[^a-zA-Z0-9]', '', names[i])
        #eliminate dups
        names_set = set(names)
        return ' '.join(names_set)

    def wall_of_alias(self):
        return ' '.join([re.sub(r'[^a-zA-Z0-9\ ]', '', alias.name.replace(',', ' ')) for alias in self.alias_set.all()])


    @property
    def entrez_str(self):
        return str(self.entrez)

    
    @property
    def ensemble_str(self):
        ensembles = []
        for g in  self.crossref_set.filter(crossrefdb__name='Ensembl'):
            ensembles.append(g.xrid)            
        return ' '.join(item for item in ensembles)



class CrossRefDB(models.Model):
    # Name of the crossref database (e.g. Entrez)
    name = models.CharField(max_length=64, unique=True, db_index=True)

    # URL of the gene page in the crossref database. The string
    # _REPL_ will be replace with the crossref gene id
    url = models.URLField()

    def natural_key(self):
        return self.name

    def __str__(self):
        return self.name

    class Meta:
        verbose_name = "Cross Reference Database"


class CrossRef(models.Model):
    crossrefdb = models.ForeignKey(CrossRefDB, on_delete=models.CASCADE)
    gene = models.ForeignKey(Gene, on_delete=models.CASCADE)
    xrid = models.CharField(max_length=32, null=False)

    def get_absolute_url(self):
        return self.crossrefdb.url.replace('_REPL_', self.xrid)

    def __str__(self):
        return self.xrid

    class Meta:
        unique_together = (("xrid", "crossrefdb"))

    def as_dict(self):
        return {'crossrefdb': self.crossrefdb.name, 'xrid': self.xrid, 'url': self.get_absolute_url()}


class Alias(models.Model):
    gene = models.ForeignKey(Gene, on_delete=models.CASCADE)
    name = models.TextField()
    organism = models.ForeignKey(Organism, on_delete=models.CASCADE)
    ordering = ['name']

    def __str__(self):
        return self.name


class GeneInterval(models.Model):
    CHROMOSOMES = (
        ('chr1', 'chromosome 1'),
        ('chr2', 'chromosome 2'),
        ('chr3', 'chromosome 3'),
        ('chr4', 'chromosome 4'),
        ('chr5', 'chromosome 5'),
        ('chr6', 'chromosome 6'),
        ('chr7', 'chromosome 7'),
        ('chr8', 'chromosome 8'),
        ('chr9', 'chromosome 9'),
        ('chr10', 'chromosome 10'),
        ('chr11', 'chromosome 11'),
        ('chr12', 'chromosome 12'),
        ('chr13', 'chromosome 13'),
        ('chr14', 'chromosome 14'),
        ('chr15', 'chromosome 15'),
        ('chr16', 'chromosome 16'),
        ('chr17', 'chromosome 17'),
        ('chr18', 'chromosome 18'),
        ('chr19', 'chromosome 19'),
        ('chr20', 'chromosome 20'),
        ('chr21', 'chromosome 21'),
        ('chr22', 'chromosome 22'),
        ('chrX', 'chromosome X'),
        ('chrY', 'chromosome Y'),

    )

    # Gene whose interval is stored
    gene = models.ForeignKey(Gene, on_delete=models.CASCADE)

    # Chromosome of gene
    chromosome = models.CharField(max_length=6, choices=CHROMOSOMES)

    # Start base pair position (inclusive) of gene variant interval
    start = models.PositiveIntegerField()

    # End base pair position (inclusive) of gene variant interval
    end = models.PositiveIntegerField()

    @property
    def just_chromosome(self):
        return self.chromosome[3:]

    @property
    def gene_entrez(self):
        return self.gene.entrez

    @property
    def difference(self):
        return int(self.end-self.start)

    class Meta:
        unique_together = ('gene', 'chromosome')

    def __str__(self):
        return '{0} ({1}:{2}) {3}'.format(self.gene, self.start, self.end, self.chromosome)
