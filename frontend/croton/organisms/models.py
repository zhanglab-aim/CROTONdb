from django.db import models
from django.template.defaultfilters import slugify

class Organism(models.Model):
    taxid       = models.PositiveIntegerField(db_index=True, unique=True, help_text="Taxonomy ID assigned by NCBI.")
    name        = models.TextField()
    slug        = models.SlugField(unique=True,help_text="Slug field, which only contains characters that are URL compatible.")

    def save(self, *args, **kwargs):
        self.slug = slugify(self.name)
        super(Organism, self).save(*args, **kwargs)

    def __str__(self):
        return self.name
