#!/bin/bash



#
# Data should be downloaded from sckidney-data-lfs
#
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
export DATADIR=$SCRIPT_DIR/../data/

function help {
cat <<EOF
Syntax: $BN [options] [DATA]
Populate scKidney from data files

DATA is one of or blank for all tables:
genes       :: Load/update genes data from NCBI entrez
index       :: Rebuild ElasticSearch index (with gene names)
Options:
-h | --help         :: this help
-d | --data-dir     :: git repository of data files (expects a /data dir)
EOF
}

function log {
    echo $(date +%F_%T) $$ $BASHPID $1
}

function die {
    echo "$BN: $@" >&2
    exit 1
}


while [ "$1" ]
do case "$1" in
    -h | --help | help ) help
        exit
	    ;;
    -d | --data-dir )
		DATADIR="$2"
	    shift 2
	    ;;
	* ) break
	    ;;
   esac
done

if [ ! -d "$DATADIR" ]
then
    die "invalid directory, not specified (DATADIR): $DATADIR"
fi

log "Importing data from $DATADIR"

cd ../croton


# Add organisms
log "Adding organism (if not added)"
python manage.py organisms_add --name="Homo sapiens" --taxid=9606


if [ -z $1 ] || [ $1 = 'genes' ]
then

    if [ ! -d "$DATADIR/genes/" ]
    then
        die "No genes directory (DATADIR/genes): $DATADIR/genes"
    fi

    #Add gene databases
    log "Adding gene databases"
    python manage.py genes_add_xrdb --name=Ensembl --URL="http://www.ensembl.org/Gene/Summary?g=_REPL_"
    python manage.py genes_add_xrdb --name=Entrez --URL="http://www.ncbi.nlm.nih.gov/gene/_REPL_"
    python manage.py genes_add_xrdb --name=HGNC --URL="https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/_REPL_"
    python manage.py genes_add_xrdb --name=HPRD --URL="http://www.hprd.org/protein/_REPL_"
    python manage.py genes_add_xrdb --name=UniProtKB --URL="http://www.uniprot.org/uniprot/_REPL_"

    #Add genes for an organism
    log "Loading geneinfo file"
    python manage.py genes_load_geneinfo --geneinfo_file=$DATADIR/genes/Homo_sapiens.gene_info --tax_id=9606 --symbol_col=2 --systematic_col=2
fi

#Gene positions from the file obtained at http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_35/gencode.v35.basic.annotation.gtf.gz
if [ -z $1 ] || [ $1 = 'positions' ]
then

    if [ ! -d "$DATADIR/genes/" ]
    then
        die "No genes directory (DATADIR/genes): $DATADIR/genes"
    fi

    GENCODE_POSITIONS=gencode.v35.basic.annotation.gtf

    if [ ! -f "$DATADIR/genes/$GENCODE_POSITIONS" ]
    then
        die "No file with gencode positions found: $DATADIR/genes/$GENCODE_POSITIONS"
    fi

    #Add genes positions 
    log "Loading gene positions file from $DATADIR/genes/$GENCODE_POSITIONS"
    python manage.py genes_load_gencode_positions --input=$DATADIR/genes/$GENCODE_POSITIONS
fi



if [ -z $1 ] || [ $1 = 'index' ]
then
    log "Rebuilding index"
    python manage.py search_index --rebuild
fi

