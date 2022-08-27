
from django.http import JsonResponse
from utils.views import HybridListView, HybridDetailView
from django.conf import settings

from .models import Gene
from .services import ElasticSearchQuerySet
from utils.hostname import resolve_hostname,timeit

class GeneDetailView(HybridDetailView):
    queryset = Gene.objects.select_related('organism')

class GeneListView(HybridListView):
    model = Gene


class GeneSearchView(HybridListView):
    '''
    GeneSearchView for accessing REST API and query for gene name.
    No querying database.
    '''

    queryset = ElasticSearchQuerySet()

    def get(self,request,*args,**kwargs):
        '''
        It queries index using term under query parameter
        :param request:
        :param args:
        :param kwargs:
        :return:
        '''
        term = request.GET.get('query')
        if not term:
            return JsonResponse([],safe=False)

        if settings.RESOLVE_HOSTNAME_FOR_SEARCH_API:
            #results = self.get_queryset().search(term,url="http://localhost:8000/search-api/genes")
            hostname = resolve_hostname(request)
            print("Hostname for gene search api: %s" % hostname)
            results = self.get_queryset().search(term,hostname=hostname)
        else:
            #Search REST API is always on the same server and on port 8000, so no need to configure hostname
            results = self.get_queryset().search(term)
        return JsonResponse(results, safe=False)


