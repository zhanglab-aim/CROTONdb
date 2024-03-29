# Generated by Django 3.2 on 2021-04-30 20:57

from django.db import migrations, models


class Migration(migrations.Migration):

    initial = True

    dependencies = [
    ]

    operations = [
        migrations.CreateModel(
            name='Organism',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('taxid', models.PositiveIntegerField(db_index=True, help_text='Taxonomy ID assigned by NCBI.', unique=True)),
                ('name', models.TextField()),
                ('slug', models.SlugField(help_text='Slug field, which only contains characters that are URL compatible.', unique=True)),
            ],
        ),
    ]
