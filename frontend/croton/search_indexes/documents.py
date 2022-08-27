from django.conf import settings

from django_elasticsearch_dsl import Document, Index, fields
from django_elasticsearch_dsl_drf.compat import KeywordField, StringField

from search_indexes.util import whitespace_lowercase_analyzer
from search_indexes.util import edge_ngram_completion_analyzer1_10, whitespace_lowercase_analyzer, \
    edge_ngram_completion_analyzer2_15, edge_ngram_completion_analyzer3_20

from genes.models import Gene

__all__ = ('GeneDocument',)

INDEX = Index(settings.ELASTICSEARCH_INDEX_NAMES[__name__])

# See Elasticsearch Indices API reference for available settings
INDEX.settings(
    number_of_shards=1,
    number_of_replicas=1,
    blocks={'read_only_allow_delete': False},
    # read_only_allow_delete=False
)

#
# To see the mapping
# http://localhost:9200/genes/_mapping
@INDEX.doc_type
class GeneDocument(Document):
#     """Gene Elasticsearch document."""
#
    entrez = StringField(attr='entrez_str',
        fields = {
            'raw':fields.KeywordField(),
            #min entrez is 1, maximum 10 digits
            'edge_ngram_completion': StringField(
                analyzer=edge_ngram_completion_analyzer1_10,
                search_analyzer=whitespace_lowercase_analyzer,
            ),
        }
    )

    standard_name = StringField(
        fields={
            'raw': KeywordField(),
            'edge_ngram_completion': StringField(
                #skipping 1 letter T gene, starting from 2 to 15 characters
                analyzer=edge_ngram_completion_analyzer2_15,
                search_analyzer=whitespace_lowercase_analyzer,
            ),
        }
    )


    systematic_name = StringField(
    fields={
        'raw': KeywordField(),
        'edge_ngram_completion': StringField(
            # skipping 1 letter T gene, starting from 2 to 15 characters
            analyzer=edge_ngram_completion_analyzer2_15,
            search_analyzer=whitespace_lowercase_analyzer,
        ),
    }
    )

    names_auto = StringField(attr='wall_of_name',
        fields={
        'edge_ngram_completion': StringField(
            analyzer=edge_ngram_completion_analyzer3_20,
            #phs-2 won't work
            #search_analyzer='standard'
            #now this works http://localhost:9999/search/genes/functional_suggest/?all_names=phs-2
            #http://localhost:9999/search/genes/functional_suggest/?all_names=rbak-LOC389458
            search_analyzer=whitespace_lowercase_analyzer,
            ),
        }
    )

    class Django(object):
        """Meta options."""
        model = Gene  # The model associate with this Document
        #don't update the index upon saving of Gene object
        ignore_signals = True
