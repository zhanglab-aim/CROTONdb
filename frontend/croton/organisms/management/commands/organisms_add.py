from django.core.management.base import BaseCommand
from organisms.models import Organism

class Command(BaseCommand):

    def add_arguments(self,parser):
        parser.add_argument('--taxid',help='Organism Taxonomy ID supplied by NCBI.')
        parser.add_argument('--name', help='Organism scientific name, e.g. "Homo sapiens"')


    help = 'Add standards from stds_file with the associations from assoc_file.'
    def handle(self, *args, **options): 
        tax_id = options['taxid']
        name = options['name']
        if tax_id and name:
        # Construct an organism object (see: organisms/models.py)
            try:
                org = Organism.objects.get(taxid=tax_id)
                org.name = name
            except Organism.DoesNotExist:
                org = Organism(taxid=tax_id, name=name)
            org.save() # Save to database specified in settings.py
            print("Organism added")
        else:
            print("Couldn't add " + str(name) + " with tax_id " + str(tax_id) + ".")
