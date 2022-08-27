from search_indexes.documents import GeneDocument
from django_elasticsearch_dsl_drf.serializers import DocumentSerializer


class GeneDocumentSerializer(DocumentSerializer):
    """Serializer for Gene document."""

    class Meta(object):
        """Meta options."""

        document = GeneDocument
        #fields need to be from GeneDocument not Gene model
        fields = (
            'entrez',
            'systematic_name',
            'standard_name',
            'names_auto',
        )
