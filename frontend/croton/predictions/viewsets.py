import os
import glob
from django.conf import settings
from django.http import Http404

from rest_framework import viewsets
from rest_framework.reverse import reverse
from rest_framework.response import Response

from genes.models import Gene, GeneInterval
from utils.hostname import resolve_hostname, timeit

from .serializers import TabixDir, TabixDirSerializer
from .serializers import GenePredictionsSerializer
from .utils import fetch_tabix_croton_predictions, convert_tabix_results, build_response_json
from .utils import CHROMOSOME_START_POSITION, DEFAULT_POSITION_RANGE, HEADER_MAP


class ChromosomePredictionsViewSet(viewsets.ViewSet):

    def list(self, request):
        dirs = []
        for dir in glob.glob('{}/*/CROTON*gz'.format(settings.TABIX_FILES_DIR)):
            chromosome = (os.path.basename(dir).split("_")[2]).split(".")[0]
            url_domain = resolve_hostname(request)
            absolute_url_path = '{}{}'.format(url_domain, reverse(
                'predictions:predictionsviewset-detail', kwargs={'pk': chromosome}))
            dirs.append(TabixDir(chromosome=chromosome,
                        path=dir, query=absolute_url_path))
        serializer = TabixDirSerializer(instance=dirs, many=True)
        # custom serializer that shows tabix paths as sell
        return Response(serializer.data)

    @timeit
    def retrieve(self, request, pk):
        chromosome = pk

        # get start/end positions for query
        start = 1
        if chromosome in CHROMOSOME_START_POSITION:
            start = CHROMOSOME_START_POSITION[chromosome]

        if request.GET.get('start'):
            start = int(request.GET.get('start'))

        start = start - 1

        end = start + DEFAULT_POSITION_RANGE
        if request.GET.get('end'):
            end = int(request.GET.get('end'))

        tabix_file = 'CROTON_varpred_{}.gz'.format(chromosome)
        tabix_file_or_url = os.path.join(os.path.join(
            settings.TABIX_FILES_DIR, chromosome), tabix_file)
        print("Reading {}".format(tabix_file_or_url))

        if not os.path.exists(tabix_file_or_url):
            print("No tabix file for chr{} found".format(chromosome))
            raise Http404("No tabix file for chr{} found".format(chromosome))

        (tabix_results, is_header) = fetch_tabix_croton_predictions(
            tabix_file_or_url,
            chrom=chromosome,
            start=start,
            end=end
        )
        header_map = {}
        # if no header included, use the default header as defined in HEADER_MAP
        if not is_header:
            header_map = HEADER_MAP

        results_len = len(tabix_results)
        print('Results length: {}'.format(results_len))

        if results_len > 0:
            response = build_response_json(tabix_results, header_map)
        else:
            response = {}
            response['results'] = {}
            response['results']['total'] = 0

        # extend response with query information
        response['query'] = {}
        response['query']['chromosome'] = chromosome
        response['query']['start'] = start
        response['query']['end'] = end

        return Response(response)


# Gene predictions per gene entrez
# e.g. /predictions/gene
# detail for gene entrez:
#   /predictions/gene/3581/
class GenePredictionsViewSet(viewsets.ReadOnlyModelViewSet):

    queryset = Gene.objects.all()
    serializer_class = GenePredictionsSerializer

    @timeit
    def retrieve(self, request, pk):
        gene_intervals = GeneInterval.objects.filter(gene__entrez=int(pk))
        # there are several cases of GeneIntervals on both chromosome X and Y but there is no predictionsn in this case for chromosome Y
        #
        gi = gene_intervals[0]
        # strip chr
        chromosome = gi.just_chromosome

        start = gi.start
        end = gi.end

        tabix_file = 'CROTON_varpred_{}.gz'.format(chromosome)
        tabix_file_or_url = os.path.join(os.path.join(
            settings.TABIX_FILES_DIR, chromosome), tabix_file)

        if not os.path.exists(tabix_file_or_url):
            print("No tabix file for chr{} found".format(chromosome))
            raise Http404("No tabix file for chr{} found".format(chromosome))

        print("\n====Querying chr{} from {} to {} for gene {}({})".format(
            chromosome, start, end, gi.gene, pk))
        # fetch results from tabix file into data frame
        (tabix_results, is_header) = fetch_tabix_croton_predictions(
            tabix_file_or_url,
            chrom=chromosome,
            start=start,
            end=end
        )
        header_map = {}
        # if no header included, use the default header as defined in HEADER_MAP
        if not is_header:
            header_map = HEADER_MAP

        response = build_response_json(tabix_results, header_map)

        # extend response with query information
        response['query'] = {}
        response['query']['chromosome'] = gi.chromosome
        response['query']['start'] = start
        response['query']['end'] = end

        return Response(response)
