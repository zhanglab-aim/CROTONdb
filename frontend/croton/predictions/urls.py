from django.urls import path, include
from rest_framework.routers import DefaultRouter
from django.views.generic.base import TemplateView
from .viewsets import (
    ChromosomePredictionsViewSet,
    GenePredictionsViewSet

)
__all__ = ('urlpatterns',)

app_name = 'predictions'

router = DefaultRouter()

router.register(r'chrom', ChromosomePredictionsViewSet, 'predictionsviewset')
router.register(r'gene', GenePredictionsViewSet, 'genepredictionsviewset')


urlpatterns = [
    path('', include(router.urls)),
    path('variant/<str>', TemplateView.as_view(
        template_name="variant.html"), name='variant')
]
