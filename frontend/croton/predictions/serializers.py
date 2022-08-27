from rest_framework import serializers
from genes.models import GeneInterval, Gene
from utils.serializers import DynamicFieldsModelSerializer

class TabixDir(object):    
    def __init__(self, **kwargs):
        for field in ('chromosome', 'path', 'query'):
            setattr(self, field, kwargs.get(field, None))

class TabixDirSerializer(serializers.Serializer):
    chromosome = serializers.CharField(max_length=3)
    path = serializers.CharField(max_length=256)    
    query = serializers.URLField()

    class Meta:
        fields = ['chromosome','path','query']


# class GeneSerializer(serializers.ModelSerializer):    
#     predictions_url = serializers.HyperlinkedIdentityField(view_name='predictions:genepredictionsviewset-detail')
#     class Meta:
#         model = Gene
#         fields = ['entrez','standard_name','gene_interval_chromosome','predictions_url']        




class GeneIntervalSerializer(serializers.ModelSerializer):
    class Meta:
        model = GeneInterval
        fields = ('start','end','chromosome')

class GenePredictionsSerializer(DynamicFieldsModelSerializer,
        serializers.HyperlinkedModelSerializer):

    url = serializers.HyperlinkedIdentityField(
        view_name = 'predictions:genepredictionsviewset-detail'
    )
    gene_intervals = GeneIntervalSerializer(many=True, read_only=True, source='get_gene_intervals')
    
    class Meta:
        model = Gene 
        fields = ('url',
                'entrez', 
                'standard_name', 
                'gene_intervals') 