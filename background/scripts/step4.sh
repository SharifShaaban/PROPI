#############################################################################
# - Need to add module and export Path for Prokka to work on server.        #
# - For each prophage annotate using the Prokka E. coli database and the    #
# faa file given (obtained from another annotation pipeline.                #
# - Copy all gbf files to separate folder.                                  #
#############################################################################
module add compilers/perl/5.18.1
export PATH=$PATH:/groups/microbial_bioinf_grp/:/usr/local/shared_bin/blast+/blast+-2.2.25/es5-64/bin/:/groups/microbial_bioinf_grp/hmmer3/:/groups/microbial_bioinf_grp/rnammer/
export PERL5LIB=/groups/microbial_bioinf_grp/perl/

raw_p_path=./temp/raw_prophages
annp_path=./temp/annotated_prophages
mkdir $annp_path &> /dev/null
mkdir ./temp/temp_gbf

rm -r $annp_path/*
rm ./temp/temp_gbf/*

for file in $raw_p_path/*.fasta
do
    filename=`echo $file| cut -d'/' -f 4`
    foldername=`echo $filename| cut -d'.' -f 1`
    /usr/local/shared_bin/prokka/prokka-1.5.2/bin/prokka --proteins ./background/panprophage.faa --genus Escherichia --species coli --usegenus --rfam --prefix $foldername --outdir $annp_path/$foldername $file --force
done

find $annp_path/*/ -name "*.gbf" -exec cp {} ./temp/temp_gbf/ \;

#############################################################################
# - Getting gene IDs from the Panprophage database for gene ontology.       #
# - DAVID gene ontology attempted in browser.                               #
#############################################################################

#cat ./temp/temp_gbf/*.gbf > ./temp/panprophage.gbf
#sed -i 's/RepID/@/g' ./temp/panprophage.gbf
#grep "@=" ./temp/panprophage.gbf | cut -d'@' -f2 | cut -d'=' -f2 | cut -d'"' -f1 > ./output/gene_IDs.txt
#sed -i 's/@/RepID/g' ./temp/panprophage.gbf

#############################################################################
# - Call python script to add colour flag to gbf file and move them to gbk  #
# folder.                                                                   #
# - Sort and unique temp file to summarise all the products falling in each #
# functional group observed.                                                #
#############################################################################

mkdir ./output/gbf &> /dev/null
rm ./output/gbf/*

anno_path=./temp/annot_groups
mkdir $annot_path &> /dev/null

python ./background/scripts/gbf_colours.py

for file in $anno_path/*.tmp
do
    filename=`echo $file| cut -d'/' -f 4`
    name=`echo $filename| cut -d'.' -f 1`
    fin_name="$anno_path/"$name".txt"
    cat $file | sort | uniq > $fin_name
    rm $file
done

#############################################################################
# - Look for IS elements in prophages using BLAST. Sort hits to have most   #
# likely and reliable IS hits.                                              #
#############################################################################

isbl_path=./temp/is_blast_res
mkdir $isbl_path &> /dev/null
rm $isbl_path/* &> /dev/null

for prophage in $raw_p_path/*.fasta
do
    filename=`echo $prophage | cut -d"/" -f4 | cut -d"." -f1`
    result_file="$isbl_path/"$filename"_is_blast.txt"
    makeblastdb -dbtype nucl -in $prophage
    blastn -db $prophage -query ./background/is_db.fasta -evalue 1e-100 -best_hit_score_edge 0.0001 -best_hit_overhang 0.25 -outfmt 6 | awk '$4>700' | sort -nrk12,12 | sort -u -k9,9 | sort -nrk12,12 | sort -u -k10,10 | sort -nk9,9 > $result_file
done

#############################################################################
# - Format the BLAST result file, to create a colour flag for IS elements   #
# containing its coordinates and adding it to the prophage GBF files.       #
#############################################################################

gbf_path=./output/gbf
mkdir $gbf_path
rm $raw_p_path/*.fasta.n* &> /dev/null
find $isbl_path/ -empty -type f -delete

for file in $isbl_path/*blast.txt
do
    sed -i "s/\t/=/g" $file
done

for file in $gbf_path/*.gbf
do
    rm ./temp_p_file.txt &> /dev/null

    g_name=`echo $file| cut -d'/' -f 4 | cut -d'_' -f 1`
    p_start=`echo $file| cut -d'/' -f 4 | cut -d'_' -f 2`
    p_end=`echo $file| cut -d'/' -f 4 | cut -d'_' -f 3`
   
    for is_line in `cat $isbl_path/$g_name*$p_end*"is_blast.txt"`
    do
        is_coo_a=`echo $is_line | cut -d'=' -f 9`
        is_coo_b=`echo $is_line | cut -d'=' -f 10`

        if [ $is_coo_a -gt $is_coo_b ]; then
            is_start=$is_coo_b
            is_end=$is_coo_a
        else
            is_start=$is_coo_a
            is_end=$is_coo_b
        fi

        printf "     IS              "$is_start".."$is_end"\n" >> temp_p_file.txt
        printf "                     /colour=255 255 0\n" >> temp_p_file.txt

    done
    
    declare -i line_number=`grep -n "ORIGIN" $file | cut -d":" -f 1`
    declare -i line_number2=$line_number-1 
    head -$line_number2 $file > ./new_file.tmp
    cat temp_p_file.txt >> ./new_file.tmp
    tail -n +$line_number $file >> ./new_file.tmp
    mv ./new_file.tmp $file
done

rm ./temp_p_file.txt &> /dev/null

for file in $isbl_path/*blast.txt
do
    sed -i "s/=/\\t/g" $file
done
