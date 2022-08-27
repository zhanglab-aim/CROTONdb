from elasticsearch_dsl import analyzer
from elasticsearch_dsl import tokenizer
from elasticsearch_dsl.analysis import token_filter


whitespace_lowercase_analyzer = analyzer(
    "whitespace_lowercase_analyzer",
    tokenizer="whitespace",
    filter=["lowercase"]
)


#EDGE NGrams 1 to 10 (e.g. for Entrez values)
edge_ngram_completion_filter1_10 = token_filter(
    'edge_ngram_completion_filter1_10',
    type="edge_ngram",
    min_gram=1,
    max_gram=10
)

#yes
edge_ngram_completion_analyzer1_10= analyzer(
    "edge_ngram_completion_analyzer1_10",
    tokenizer="whitespace",
    filter=["lowercase", edge_ngram_completion_filter1_10]
)


#EDGE NGrams 2 to 15 (e.g. for gene symbols)
edge_ngram_completion_filter2_15 = token_filter(
    'edge_ngram_completion_filter2_15',
    type="edge_ngram",
    min_gram=2,
    max_gram=15
)

edge_ngram_completion_analyzer2_15= analyzer(
    "edge_ngram_completion_analyzer2_15",
    tokenizer="whitespace",
    filter=["lowercase", edge_ngram_completion_filter2_15]
)


#min_gram 3 to 20
edge_ngram_completion_filter3_20 = token_filter(
    'edge_ngram_completion_filter3_20',
    type="edge_ngram",
    min_gram=1,
    max_gram=20

)

edge_ngram_completion_analyzer3_20 = analyzer(
    "edge_ngram_completion_analyzer3_20",

    #tokenize by whitespace, so BTO:00045 will be preserved
    tokenizer="whitespace",
    filter=["lowercase", edge_ngram_completion_filter3_20]
)
