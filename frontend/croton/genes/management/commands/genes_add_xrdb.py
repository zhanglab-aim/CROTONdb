import logging
logger = logging.getLogger(__name__)

from django.core.management.base import BaseCommand
from genes.models import CrossRefDB

class Command(BaseCommand):
    def add_arguments(self,parser):
        parser.add_argument('--name',help='Database name')
        parser.add_argument('--URL', help='Database url')

    help = 'Add a crossreference database'
    def handle(self, *args, **options):
        name = options['name']
        url = options['URL']
        if name and url:
            xrdb, created = CrossRefDB.objects.get_or_create(name=name)
            xrdb.url = url
            xrdb.save()
            print("Created")
        else:
            print("Please provide --name and --URL")
