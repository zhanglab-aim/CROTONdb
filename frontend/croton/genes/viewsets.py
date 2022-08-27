from django.shortcuts import get_object_or_404

from rest_framework import viewsets
from rest_framework import mixins
from rest_framework.response import Response
from rest_framework.decorators import detail_route

import logging
logger = logging.getLogger(__name__)

from .models import Gene 
from .serializers import GeneSerializer 

from predictions.models import Prediction
from predictions.serializers import PredictionSerializer
from terms.models import Database, Annotation
from terms.serializers import TermSerializer, AnnotationSerializer

class GeneViewSet(mixins.RetrieveModelMixin, 
        mixins.ListModelMixin, 
        viewsets.GenericViewSet):
    queryset = Gene.objects.all()
    serializer_class = GeneSerializer

    @detail_route()
    def predictions(self, request, pk=None):
        """
        Returns a list of terms a gene is predicted to be associated with 

        URL parameters:
            version: verison identifier of predictions 
            database: database of terms to filter predictions 
            score_cutoff: minimum score of predictions to return
        """
        predictions = Prediction.objects.filter(gene=pk).order_by('-score')

        version = self.request.query_params.get('version', None)
        if version:
            predictions = predictions.filter(version=version)

        database = self.request.query_params.get('database', None)
        if database:
            database = get_object_or_404(Database, slug=database)
            predictions = predictions.filter(term__database=database)

        score_cutoff = self.request.query_params.get('score_cutoff',None)
        if score_cutoff:
            predictions = predictions.filter(score__gt=score_cutoff)

        serializer = PredictionSerializer(predictions, 
                fields=('score','term'), 
                many=True, 
                context={'request':request})

        return Response(serializer.data)


    @detail_route()
    def annotations(self, request, pk=None):
        """
        Returns a list of terms a gene is annotated to 

        URL parameters:
            database: database of terms to filter predictions 
            max_term_size: filter terms that have more than max_term_size
                            gene annotations
        """

        gene = get_object_or_404(self.queryset, pk=pk)

        annotations = Annotation.objects.filter(gene=gene)

        database = self.request.query_params.get('database', None)
        if database:
            database = get_object_or_404(Database, slug=database)
            annotations = annotations.filter(term__database=database)

        max_term_size = self.request.query_params.get('max_term_size', None)
        if max_term_size and max_term_size.isdigit():
            annotations = [annot for annot in annotations if annot.term.get_total_annotations() < int(max_term_size)]

        serializer = AnnotationSerializer(annotations, 
                many=True, 
                context={'request':request})

        return Response(serializer.data)
