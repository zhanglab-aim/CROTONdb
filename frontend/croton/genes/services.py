import requests
from django.urls import reverse

def get_unique_response(url,query_str,term):
    '''
    :param url: URL to query
    :param query_str: query parameter
    :param token:  query token (search word)
    :param type_id_set: sets of already found ids
    :param results: list of results
    :return:
    '''

    results = []
    type_id_set = set()
    print(url)
    print(term)
    params = {query_str: term, 'format': 'json'}
    try:
        r = requests.get(url, params=params)
    except requests.exceptions.RequestException as e:
        print(e)
        return ['Connection error']
    response_json = r.json()
    
    print(len(response_json["results"]))
    for entry in response_json["results"]:
        new_dict = {}
        new_dict["type"] = "gene"
        new_dict["id"] = entry["entrez"]
        new_dict["symbol"] = entry["standard_name"]
        new_dict["name"] = entry["systematic_name"]
        new_dict["aliases"] = entry["names_auto"]

        #do not add dups
        if not new_dict["id"] in type_id_set:
            type_id_set.add(new_dict["id"])
            results.append(new_dict)

    return results


class ElasticSearchQuerySet(object):
    '''Class to represent searching ES '''

    def __init__(self,host=None):
        self.host = "http://localhost:8000"
        if host:
            self.host = host

        self.query_str = "search_multi_match"

    def search(self,term,url=None,hostname=None,query_str=None):


        if hostname:
            self.host = hostname

        self.url = self.host + reverse('search:genedocument-list')

        if url:
            self.url = url

        if query_str:
            self.query_str=query_str

        print(self.query_str)
        results_list = []

        results_list = get_unique_response(self.url, self.query_str, term)

        return results_list
