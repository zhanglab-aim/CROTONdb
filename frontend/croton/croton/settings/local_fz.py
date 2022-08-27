'''
Local settings
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

# docker
DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.mysql',
        'NAME': 'croton',
        'USER': 'root',
        'PASSWORD': 'mysql123',
        'HOST': '127.0.0.1',

    }
}
# print(DATABASES)

# SECRET CONFIGURATION
# ------------------------------------------------------------------------------
# See: https://docs.djangoproject.com/en/dev/ref/settings/#secret-key
# Note: This key only used for development and testing.
SECRET_KEY = 'n###@3(uzdf#!2@+&amp;htdaafsf==@vk2Mit$kyNo*odj-'


# This is for ES6 docker run -p 9200:9200 -p 9300:9300 -e "discovery.type=single-node" docker.elastic.co/elasticsearch/elasticsearch:6.2.4
STATICFILES_DIRS = [os.path.join(BASE_DIR, "static")]
ELASTICSEARCH_DSL = {
    'default': {
        'hosts': 'localhost:9200'
    },
}

# Name of the Elasticsearch index
ELASTICSEARCH_INDEX_NAMES = {
    'search_indexes.documents': 'test_gene',
}

MEDIA_URL = '/media/'
MEDIA_ROOT = "/home/atadych/dev/workspace/python3/croton_project/croton/media/"
TABIX_FILES_DIR = '/Users/alicja/dev/data/croton/tabix'


TEST_ECDF_FILE = '/Users/alicja/dev/data/croton/test_gw_data.pkl'

# REST_FRAMEWORK = {
#     'DEFAULT_PAGINATION_CLASS': 'predictions.paginator.CustomResultsSetPagination',
# }
