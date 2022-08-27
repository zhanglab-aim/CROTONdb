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

#mysql -P 4406 --protocol=tcp -u root -p
#docker
DATABASES = {
    'default' : {
    'ENGINE': 'django.db.backends.mysql',
    'NAME': 'sckidney',
    #mysql settings the same as in mysql/secrets.env
    'USER': 'root',
    'PASSWORD': 'mysql123',

    #docker-dev-mod
    #'PASSWORD': 'mysql',
    #port as in docker-compose.yml
    #'PORT': 4406,
    #docker-dev-mod
    #'PORT': 43333,

    }
}



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
BNSERVER_PORT=1234

# Prior probability for bayesian
PRIOR = 0.1

# Public gene sets pickle
#PGENES = PROJECT_ROOT + '/pickles/public_genes.p'


STATICFILES_DIRS = [os.path.join(BASE_DIR,"static")]

ELASTICSEARCH_DSL={
    'default': {
    #docker-dev-mod
    'hosts': 'localhost:49200'
    #    'hosts': 'localhost:9200'
    },
}

# Name of the Elasticsearch index
ELASTICSEARCH_INDEX_NAMES = {
    'search_indexes.documents': 'test_gene',
}
MEDIA_ROOT="/Users/alicja/dev/workspace/scKidney/media/"