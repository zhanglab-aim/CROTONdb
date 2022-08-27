def resolve_hostname(request):
    '''
    Method resolves hostname and protocol to build a correct URL
    :param request:
    :return:
    '''
    if request.is_secure():
        protocol = 'https'
    else:
        protocol = 'http'
    hostname = protocol + "://" + request.get_host()
    parts = hostname.split(":")
    if len(parts) < 3:
        hostname = hostname +":"+request.get_port()
    return hostname

import time
def timeit(method):
    def timed(*args, **kw):
        ts = time.time()
        result = method(*args, **kw)
        te = time.time()
        if 'log_time' in kw:
            name = kw.get('log_name', method.__name__.upper())
            kw['log_time'][name] = int((te - ts) * 1000)
        else:
            print(('%r  %2.2f ms') % \
                  (method.__name__, (te - ts) * 1000))
        return result
    return timed