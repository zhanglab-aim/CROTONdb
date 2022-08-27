'''
Local settings for initializing services with docker
'''


from .base import *
import os

ADMINS = (
     ('Alicja Tadych', 'atadych@princeton.edu')
)

MANAGERS = ADMINS

# DEBUG
# ------------------------------------------------------------------------------
DEBUG = True
TEMPLATE_DEBUG = DEBUG


#local mysql
# #mysql -u root -p
DATABASES = {
    'default' : {
    'ENGINE': 'django.db.backends.mysql',
    #'ENGINE': 'mysql.connector.django',
    'NAME': 'sckidney',
    #mysql settings the same as in mysql/secrets.env
    'USER': 'root',
    'PASSWORD': 'mysql123',
    'HOST': '127.0.0.1',
    #port as in docker-compose.yml
    #'PORT': 3333,

    }
}
#docker
#mysql -P 43333 --protocol=tcp -u root -p
# DATABASES = {
#     'default' : {
#     'ENGINE': 'django.db.backends.mysql',
#     #'ENGINE': 'mysql.connector.django',
#     'NAME': 'sckidney',
#     #mysql settings the same as in mysql/secrets.env
#     'USER': 'root',
#     'PASSWORD': 'mysql',
#     #port as in docker-compose.yml
#     #'HOST': 'mysql',
#     'PORT': 43333,
#
#     }
# }


# SECRET CONFIGURATION
# ------------------------------------------------------------------------------
# See: https://docs.djangoproject.com/en/dev/ref/settings/#secret-key
# Note: This key only used for development and testing.
SECRET_KEY = 'ns$jd@asfd&^zdf#!2@+&amp;htvk2Mit$kyNo*odj-'


#
# BNServer settings
#
BNSERVER_HOST="localhost"
PRIOR=0.1
#docker setup
BNSERVER_PORT=41234
#BNSERVER_PORT=9999
#ssh -L 9999:localhost:9999 atadych@crucio

# Prior probability for bayesian
PRIOR = 0.1

# Public gene sets pickle
#PGENES = PROJECT_ROOT + '/pickles/public_genes.p'


STATICFILES_DIRS = [os.path.join(BASE_DIR,"static")]

ELASTICSEARCH_DSL={
    'default': {
        'hosts': 'localhost:49200'
    },
}

# Name of the Elasticsearch index
ELASTICSEARCH_INDEX_NAMES = {
    'search_indexes.documents': 'test_gene',
}

#Settings which Gene Expression version should be taken (table expression_version id)
GENE_EXPRESSION_VERSION_ID=2
MEDIA_ROOT="/Users/alicja/dev/workspace/scKidney/media/"