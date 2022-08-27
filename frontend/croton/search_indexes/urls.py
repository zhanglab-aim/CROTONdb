from django.urls import path,include
from rest_framework.routers import DefaultRouter

from .viewsets import (
    GeneDocumentViewSet
)
__all__ = ('urlpatterns',)

app_name = 'search'

router = DefaultRouter()

genes = router.register(
    r'genes',
    GeneDocumentViewSet,
    basename='genedocument'
)


urlpatterns = [
    path('', include(router.urls)),
]
