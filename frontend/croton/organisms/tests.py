from django.test import TestCase
from django.template.defaultfilters import slugify
from django.core import management

from organisms.models import Organism

class OrganismTest(TestCase):
    def setUp(self):
        Organism.objects.create(taxid=9606, name="Homo sapiens")

    def tearDown(self):
        Organism.objects.all().delete()

    def test_slugs(self):
        """
        Tests that save creates the organism's slug.
        """
        human = Organism.objects.get(taxid=9606)
        self.assertEqual(human.slug, slugify("Homo sapiens"))

    def test_management(self):
        """
        Tests the management command that creates organisms.
        """
        management.call_command('organisms_add', taxid=6238, name="Caenorhabditis elegans")
