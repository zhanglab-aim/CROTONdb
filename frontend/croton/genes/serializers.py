from rest_framework import serializers

from .models import Gene, CrossRef, CrossRefDB, GeneInterval
from utils.serializers import DynamicFieldsModelSerializer

class CrossRefDBSerializer(serializers.ModelSerializer):
    class Meta:
        model = CrossRefDB
        fields = ('name', 'url')

class CrossRefSerializer(serializers.ModelSerializer):
    crossrefdb = CrossRefDBSerializer(read_only=True)
    class Meta:
        model = CrossRef
        fields = ('xrid','crossrefdb')


class GeneSerializer(DynamicFieldsModelSerializer,
        serializers.HyperlinkedModelSerializer):

    url = serializers.HyperlinkedIdentityField(
        view_name = 'api:gene-detail'
    )
    aliases = serializers.ReadOnlyField(source='get_aliases')
    xrefs = CrossRefSerializer(many=True, read_only=True, source='get_xrefs')
    
    class Meta:
        model = Gene 
        fields = ('url',
                'entrez', 
                'systematic_name',
                'standard_name', 
                'description', 
                'xrefs',
                'aliases') 