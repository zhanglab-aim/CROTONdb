from django.urls import path, include
from django.contrib import admin
from django.views.generic.base import TemplateView
from django.conf import settings
from django.conf.urls.static import static

urlpatterns = [
    path('', TemplateView.as_view(template_name="home.html"), name='home'),
    path('about', TemplateView.as_view(template_name="about.html"), name='about'),
    path('help', TemplateView.as_view(template_name="help.html"), name='help'),
    path('api', TemplateView.as_view(template_name="api.html"), name='api'),
    path('search-api/', include('search_indexes.urls', namespace='search')),
    path('genes/', include('genes.urls', namespace='genes')),
    path('predictions/', include('predictions.urls', namespace='predictions')),
    path('admin/', admin.site.urls),

]

urlpatterns += static(settings.STATIC_URL, document_root=settings.STATIC_ROOT)
urlpatterns += static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)
