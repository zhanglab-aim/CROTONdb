#!/usr/bin/env bash
#SBATCH --time=2-5
#SBATCH --partition=ccb
#VL 080621
#for i in `seq 1 12`; do echo $i; sbatch create_files.sh $i; done 
#sbatch create_files.sh X
echo `date`
echo begin $0
# $1 represents inputted chromosome, use start--> :' end --> ' to comment out section

cd /mnt/home/vli/ceph/croton

# Run make_beds.py
python src/variant/make_beds.py --chrom $1

# Run make_tsvs.py
if [[ $? = 0 ]]; 
    then 
        echo "Finished making beds!"
       python src/variant/make_tsvs.py --chrom $1
else
    echo "Error: make_beds did not run successfully" 
    exit 1
fi

# Run make_tabix.py
if [[ $? = 0 ]]; 
    then echo "Finished making tsvs!" 
    python src/variant/make_tabix.py --chrom $1
else 
    echo "Error: make_tsvs did not run successfully" 
    exit 1
fi

# Bgzip created vcf file
if [[ $? = 0 ]]; 
    then echo "Finished making vcf file!" 
    cd /mnt/home/vli/ceph/croton/datavl/variant/tabix
    bgzip -c "$1"/CROTON_varpred_"$1".tsv > "$1"/CROTON_varpred_"$1".gz
else 
    echo "Error: bgzip did not run successfully"
    exit 1
fi

# Tabix created bgziped '.gz' file --> 'gz.tbi'
if [[ $? = 0 ]]; 
    then echo "Finished zipping and now making tabix...!!!" 
    tabix -p vcf "$1"/CROTON_varpred_"$1".gz
else 
    echo "Error making tabix file was unsuccessful" 
    exit 1
fi
