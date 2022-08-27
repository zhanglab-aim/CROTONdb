"""Production settings and globals."""

from __future__ import absolute_import

from os import environ

from .base import *

# Normally you should not import ANYTHING from Django directly
# into your settings, but ImproperlyConfigured is an exception.
from django.core.exceptions import ImproperlyConfigured


def get_env_setting(setting):
    """ Get the environment setting or return exception """
    try:
        return environ.get(setting)
    except KeyError:
        error_msg = "Set the %s env variable" % setting
        raise ImproperlyConfigured(error_msg)

########## HOST CONFIGURATION
# See: https://docs.djangoproject.com/en/1.5/releases/1.5/#allowed-hosts-required-in-production
ALLOWED_HOSTS = ['croton.princeton.edu','croton','localhost']
########## END HOST CONFIGURATION

STATICFILES_DIRS = [os.path.join(BASE_DIR,"static")]
print(STATICFILES_DIRS)
STATIC_ROOT="/Genomics/local/croton/www/static"
MEDIA_ROOT="/Genomics/local/croton/www/media"

DATABASES = {
    "default": {
        "ENGINE": "django.db.backends.mysql",
        "NAME": "croton",
        "USER": "croton",
        "PASSWORD":get_env_setting('DB_PASSWORD'),
        "PORT": "3306",
    },
}

#DATABASES = {
#    'default': {
#        'ENGINE': 'django.db.backends.sqlite3',
#        'NAME':  '/Genomics/local/croton/www/croton/croton.sqlite3',
#    }
#}



########## SECRET CONFIGURATION
# See: https://docs.djangoproject.com/en/dev/ref/settings/#secret-key
SECRET_KEY ="asf$AFBg%9=q#xbqreedekg8^a!%xa*&7zduf0^pyss$6b(vm!9h#s"


ELASTICSEARCH_DSL={
    'default': {
        'hosts': 'gen-501r.princeton.edu:9200'
    },
}

# Name of the Elasticsearch index
ELASTICSEARCH_INDEX_NAMES = {
    'search_indexes.documents': 'croton_gene',
}
#Needs to resolve hostname, so can it access it using https://croton.princeton.edu:443
RESOLVE_HOSTNAME_FOR_SEARCH_API=True


DEBUG=False





# Google analytics
GOOGLE_ANALYTICS_KEY = 'UA-121663754-1?????'
TABIX_FILES_DIR='/express/croton/tabix'
