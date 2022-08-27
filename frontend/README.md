# croton_project #

### To use this project follow these steps: ###

* Check the project requirements (python3, etc)
* Clone a project from github
* Create working environment
* Edit settings
* Create database
* Run initial migrations 
* Populate databases

### Requirements ###
* Install `pipenv` (https://pipenv.readthedocs.io/en/latest/) and `python3`
* Install database server (MySQL)
* Run ElasticSearch (via docker).  

        $ docker run -p 9200:9200 -p 9300:9300 -e "discovery.type=single-node" docker.elastic.co/elasticsearch/elasticsearch:6.2.4




## Clone a project ##

Clone a project from the bitbucket:

    $ git clone https://github.com/FunctionLab/croton_project.git

If you experience problems with authentication, make sure SSH keys are defined.

    

## Working Environment (using pipenv)##

    $ cd croton_project/croton
If you want to have virtualenv in the project directory, export these variables
    
    $ export PIPENV_VENV_IN_PROJECT=TRUE
    $ export PIPENV_MAX_DEPTH=1
    
Now install Django and all requirements (from Pipfile.lock)

    $ pipenv install    

    
To use the environment      

    $ pipenv shell


## Using docker to initialize ElasticSearch and MySQL

Elastic Search v6   

        $ docker run -p 9200:9200 -p 9300:9300 -e "discovery.type=single-node" docker.elastic.co/elasticsearch/elasticsearch:6.2.4

## Installing all the services on the same localhost

### Database
    
#### Edit settings ##

Depending on where your setup is, edit `local.py` in croton/settings/. Make sure that DJANGO_SETTINGS_MODULE points to the desired configuration:

    $ export DJANGO_SETTINGS_MODULE=croton.settings.local

#### Create a database ##

Create a database corresponding to the name specified in the settings file under 'DATABASES' with the necessary username and password.

#### Run initial migrations ##

Run the migrations that will install necessary tables. 
All the migrations are in the repository, so you just need to apply them.

Run command:
 
    $ python manage.py migrate


#### Run tests ##

In predictions (to test tabix file):

    $ python manage.py test predictions
    
All tests:  

    $ python manage.py test


#### Populate database and create search index ##
Change DATADIR in setup.sh to point to a directory where genes file is `Homo_sapiens.gene_info` (from ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz)
Gene positions file e.g. `gencode.v35.basic.annotation.gtf` should be downloaded from  http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_35/

Run command:

    $ cd scripts
    $ setup.sh


#### Run local webserver

    $ python manage.py runserver

and check in the webserver `localhost:8000`
