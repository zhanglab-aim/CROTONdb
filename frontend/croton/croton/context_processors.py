from genes.models import GeneInterval
from django.conf import settings


class Chromosome:
    pass


def chromosomes(request):
    print("chromosomes")
    lst = []
    for g in GeneInterval.CHROMOSOMES:
        chr = Chromosome()
        chr.id = g[0][3:]
        chr.label = g[1]
        lst .append(chr)
    return {'chromosomes': lst}


def ga_key(request):
    if not hasattr(settings, 'GOOGLE_ANALYTICS_KEY'):
        print("====GOOGLE ANALYTICS not set")
        return {}
    print("===GOOGLE ANALYTICS is set")
    # add the google analytics key to the context
    return {'GOOGLE_ANALYTICS_KEY': settings.GOOGLE_ANALYTICS_KEY}
