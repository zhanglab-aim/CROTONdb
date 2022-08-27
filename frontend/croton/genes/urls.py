from django.urls import path
from . import views

app_name = 'genes'

urlpatterns = [
    path('', views.GeneListView.as_view(), name="list"),
    path('<int:pk>/', views.GeneDetailView.as_view(), name="detail"),
    #eg. http://localhost:8000/genes/search/?query=brca
    path('search/', views.GeneSearchView.as_view(), name="search"),
]

