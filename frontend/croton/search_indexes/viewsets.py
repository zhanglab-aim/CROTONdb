from rest_framework import permissions
from rest_framework.decorators import permission_classes

from django_elasticsearch_dsl_drf.filter_backends import (
    FilteringFilterBackend,
    DefaultOrderingFilterBackend,
    OrderingFilterBackend,
    MultiMatchSearchFilterBackend
)
from django_elasticsearch_dsl_drf.viewsets import DocumentViewSet

from search_indexes.documents import GeneDocument
from search_indexes.serializers import GeneDocumentSerializer
from django_elasticsearch_dsl_drf.pagination import LimitOffsetPagination


@permission_classes((permissions.AllowAny,))
class GeneDocumentViewSet(DocumentViewSet):
    """The GeneDocumentViewSet view."""

    document = GeneDocument
    serializer_class = GeneDocumentSerializer

    pagination_class = LimitOffsetPagination
    lookup_field = 'entrez'
    filter_backends = [
        FilteringFilterBackend,
        OrderingFilterBackend,
        DefaultOrderingFilterBackend,
        MultiMatchSearchFilterBackend
    ]

    # Define filtering fields
    # Necessary if you want to search by field name
    #http://localhost:8000/search/genes/?entrez=123
    #http://localhost:8000/search/genes/?standard_name=ERV3-1
    filter_fields = {
        'entrez': None,
        'standard_name': 'standard_name.raw',
    }
    # Define ordering fields
    ordering_fields = {
        'entrez':None
        #'standard_name': None,
    }
    # Specify default ordering
    ordering = ('_score','standard_name', )

    #http://localhost:8000/search/genes/?search_multi_match=123
    #http://localhost:8000/search/genes/?search_multi_match=CCR1
    multi_match_search_fields = {
        'entrez': {'boost': 4},
        'entrez.edge_ngram_completion': {'boost': 3},
        'standard_name': {'boost': 3},
        'standard_name.edge_ngram_completion': {'boost': 2},
        'names_auto.edge_ngram_completion': None,
    }