'''
Local settings

- Run in Debug mode
- Add Django Debug Toolbar
- Add django-extensions as app
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

# DATABASE
# ------------------------------------------------------------------------------
#default
# DATABASES = {
#     'default' : {
#     'ENGINE': 'django.db.backends.mysql',
#     'NAME': 'sckidney',
#     'USER': 'root',
#     'PASSWORD': '',
#     'HOST': '',
#     }
# }
#docker
DATABASES = {
    'default' : {
    'ENGINE': 'django.db.backends.mysql',
    'NAME': 'sckidney',
        "USER": "sckidney",
        "PASSWORD": "silly",
        "HOST": "/var/lib/mysql/mysql.sock",
        "PORT": "3306",
    },
}


#print(DATABASES)

# SECRET CONFIGURATION
# ------------------------------------------------------------------------------
# See: https://docs.djangoproject.com/en/dev/ref/settings/#secret-key
# Note: This key only used for development and testing.
SECRET_KEY = 'ns$jd@3(uzdf#!2@+&amp;htvk2Mit$kyNo*odj-'


# django-debug-toolbar
# ------------------------------------------------------------------------------
#IDDLEWARE_CLASSES = ('debug_toolbar.middleware.DebugToolbarMiddleware',) + MIDDLEWARE_CLASSES
#INSTALLED_APPS += ('debug_toolbar', )
#INTERNAL_IPS = ('127.0.0.1', '10.0.2.2',)


# MIDDLEWARE_CLASSES =  ('corsheaders.middleware.CorsMiddleware',
#     'django.middleware.common.CommonMiddleware') + MIDDLEWARE_CLASSES
# MIDDLEWARE_CLASSES = ('django.middleware.cache.UpdateCacheMiddleware',
#     'django.middleware.common.CommonMiddleware',
#     'django.middleware.cache.FetchFromCacheMiddleware') + MIDDLEWARE_CLASSES
#
# CORS_ORIGIN_ALLOW_ALL = True



# Settings for vegas module. URL specification is necessary to get clickable
# links in a POST response. (GET calls correctly prepend URLs for uploaded files
# even without this but POST calls don't.)
#MEDIA_ROOT = os.path.join(PROJECT_ROOT, 'media/')
#MEDIA_URL  = 'http://scda002.scdanet.org:8000/media/'

# Maximum upload size of the Vegas file (10MB)
MAX_UPLOAD_SIZE = "10485760"

#
# BNServer settings
#
#ssh -L 9999:localhost:9999 crucio
#giant v1
#BNSERVER_PORT=1234
#giant v2
BNSERVER_PORT=9999
BNSERVER_HOST="portus.princeton.edu"
PRIOR=0.1



# Prior probability for bayesian
PRIOR = 0.1

# Public gene sets pickle
#PGENES = PROJECT_ROOT + '/pickles/public_genes.p'


STATICFILES_DIRS = [os.path.join(BASE_DIR,"static")]
ELASTICSEARCH_DSL={
    'default': {
        'hosts': 'gen-501r.princeton.edu:9200'
    },
}

# Name of the Elasticsearch index
ELASTICSEARCH_INDEX_NAMES = {
    'search_indexes.documents': 'gene',
}

MEDIA_URL='/media/'
MEDIA_ROOT="/Genomics/local/sckidney/www/media"
MEDIA_ROOT="/var/www/sckidney/www/media"
STATIC_ROOT="/Genomics/local/sckidney/www/static"
STATIC_ROOT="/var/www/sckidney/www/static"

ALLOWED_HOSTS=["sckidney.princeton.edu"]
#print(TEMPLATES)
CELL_TRACJECTORY_SVG_IMAGE='/var/www/sckidney/www/media/viz/cell_trajectory_bare_image_V2.svg'
ANATOMICAL_MAP_SVG_IMAGE='/var/www/sckidney/www/media/viz/anatomical_map_bare_image.svg'
#ANATOMICAL_MAP_SVG_IMAGE='/Genomics/local/sckidney/www/media/viz/anatomical_map_bare_image.svg'
#CELL_TRACJECTORY_SVG_IMAGE='/Genomics/local/sckidney/www/media/viz/cell_trajectory_bare_image.svg'

GENE_EXPRESSION_VERSION_ID=2

#Pretubular aggregate 2, find integration id for this context
#select i.id as integration_id, i.context_id, c.title from integrations_integration i join contexts_context c on i.context_id = c.id where title like 'Pretu%2';i
DEF_NET_INTEG_CELL_TYPE_ID=172
DEF_NET_INTEG_CELL_LINEAGE_ID=166
DEF_NET_INTEG_GLOBAL_ID=189


#Needs to resolve hostname, so can it access it using https://sckidney.princeton.edu:443
RESOLVE_HOSTNAME_FOR_SEARCH_API=True


GOOGLE_ANALYTICS_KEY='UA-165269676-1'
