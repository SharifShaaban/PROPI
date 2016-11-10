module add compilers/python/2.7.9_shared
module add apps/gcc/R/3.2.2
easyfig_path=/groups2/gally_grp/fsa_project/sharif/tools/Easyfig_CL_2.1.py

#############################################################################
# - Calculates the average prophage length, and the average number of genes #
# of a specific functional group using their colour flags. Calculates the   #
# proportion of those genes for each cluster, but multiplies the proportion #
# with the maximum prophage length of all prophages to have numbers that can#
# be worked with.                                                           #
# scale used for decimal significant number, bc for mathematical calculation#
#############################################################################

group_path=./temp/groupings

rm ./temp/easyfig_rc_log.txt &> /dev/null
rm ./temp/*length_log.txt &> /dev/null
mkdir ./temp/length_logs/ &> /dev/null
rm ./temp/length_logs/* &> /dev/null
mkdir ./output/easyfig_res &> /dev/null
rm ./output/easyfig_res/* &> /dev/null

declare -i max_length=`ls -l ./output/gbf/*.gbf | cut -d'/' -f4 | cut -d'_' -f4 | cut -d'.' -f1 | sort -n | tail -1`

for folder in `cat ./temp/group_lst/*_4.5.txt`
do
    name=`echo $folder| cut -d'/' -f 1`
    group_name=$name
    declare -i counter=0
    declare -i full_length=0
    for file in ./temp/groupings/$name/*
    do
        prophagename=`echo $file | cut -d'/' -f 5 | cut -d'.' -f 1`
        length=`echo ${prophagename##*_}`
        declare -i int_length=$length
        full_length=$full_length+$int_length
        counter=$counter+1
    done   

    declare -i average_length=$full_length/$counter

    declare -i all_LYS_genes=`grep "/colour=32 160 64" $group_path/$name/* | wc -l`
    average_LYS_genes=`echo "scale=2; $all_LYS_genes / $counter" | bc -l`
    LYS_proportion=`echo "scale=2; $average_LYS_genes * $max_length / $average_length" | bc -l`
    printf "\n"$name"\t"$average_length"\t"$average_LYS_genes"\t"$LYS_proportion"\n" >> LYS_length_log.txt

    declare -i all_REG_genes=`grep "/colour=33 96 192" $group_path/$name/* | wc -l`
    average_REG_genes=`echo "scale=2; $all_REG_genes / $counter" | bc -l`
    REG_proportion=`echo "scale=2; $average_REG_genes * $max_length / $average_length" | bc -l`
    printf "\n"$name"\t"$average_length"\t"$average_REG_genes"\t"$REG_proportion"\n" >> REG_length_log.txt

    declare -i all_EFF_genes=`grep "/colour=224 128 128" $group_path/$name/* | wc -l`
    average_EFF_genes=`echo "scale=2; $all_EFF_genes / $counter" | bc -l`
    EFF_proportion=`echo "scale=2; $average_EFF_genes * $max_length / $average_length" | bc -l`
    printf "\n"$name"\t"$average_length"\t"$average_EFF_genes"\t"$EFF_proportion"\n" >> EFF_length_log.txt

    declare -i all_REC_genes=`grep "/colour=166 202 240" $group_path/$name/* | wc -l`
    average_REC_genes=`echo "scale=2; $all_REC_genes / $counter" | bc -l`
    REC_proportion=`echo "scale=2; $average_REC_genes * $max_length / $average_length" | bc -l`
    printf "\n"$name"\t"$average_length"\t"$average_REC_genes"\t"$REC_proportion"\n" >> REC_length_log.txt

    declare -i all_STR_genes=`grep "/colour=160 224 128" $group_path/$name/* | wc -l`
    average_STR_genes=`echo "scale=2; $all_STR_genes / $counter" | bc -l`
    STR_proportion=`echo "scale=2; $average_STR_genes * $max_length / $average_length" | bc -l`
    printf "\n"$name"\t"$average_length"\t"$average_STR_genes"\t"$STR_proportion"\n" >> STR_length_log.txt

    declare -i all_TRN_genes=`grep "/colour=96 32 128" $group_path/$name/* | wc -l`
    average_TRN_genes=`echo "scale=2; $all_TRN_genes / $counter" | bc -l`
    TRN_proportion=`echo "scale=2; $average_TRN_genes * $max_length / $average_length" | bc -l`
    printf "\n"$name"\t"$average_length"\t"$average_TRN_genes"\t"$TRN_proportion"\n" >> TRN_length_log.txt

    declare -i all_MET_genes=`grep "/colour=192 160 192" $group_path/$name/* | wc -l`
    average_MET_genes=`echo "scale=2; $all_MET_genes / $counter" | bc -l`
    MET_proportion=`echo "scale=2; $average_MET_genes * $max_length / $average_length" | bc -l`
    printf "\n"$name"\t"$average_length"\t"$average_MET_genes"\t"$MET_proportion"\n" >> MET_length_log.txt

    declare -i all_OTH_genes=`grep "/colour=192 192 192" $group_path/$name/* | wc -l`
    average_OTH_genes=`echo "scale=2; $all_OTH_genes / $counter" | bc -l`
    OTH_proportion=`echo "scale=2; $average_OTH_genes * $max_length / $average_length" | bc -l`
    printf "\n"$name"\t"$average_length"\t"$average_OTH_genes"\t"$OTH_proportion"\n" >> OTH_length_log.txt

    declare -i all_IS=`grep "/colour=255 255 0" $group_path/$name/* | wc -l`
    average_IS=`echo "scale=2; $all_IS / $counter" | bc -l`
    IS_proportion=`echo "scale=2; $average_IS * $max_length / $average_length" | bc -l`
    printf "\n"$name"\t"$average_length"\t"$average_IS"\t"$IS_proportion"\n" >> IS_length_log.txt

    printf "\n"$name"\t"$counter"\t"$average_length"\n" >> summary_length_log.txt 

#############################################################################
# - If loops created to sort between singleton clusters and others.         #
# - In the case of singleton clusters their sole member is added to the log #
# file with its number of genes of a functional group and its length, as    #
# opposed to the averages and proportions previously recorded for each      #
# cluster (this is at an individual prophage level).                        #
# - This is also done for clusters with more than one member but an extra   #
# step is included to determine whether some prophages require to be reverse#
# complemented (using the EMBOSS revseq tool) so that all prophages in the  #
# cluster are orientated in the same manner for the Easyfig figures which   #
# are created further on (Reverse complemented prophages recorded in        #
# "easyfig_rc_log.txt".                                                     #
# - At the end of the loop command line Easyfig is used to make cluster     #
# prophage alignments.                                                      #
#############################################################################

    feature_files=""
    if [ $(ls -l ./temp/groupings/$name/* | wc -l) -gt 1 ];then
        for file in ./temp/groupings/$name/*
        do
            filename=`echo $file | cut -d'/' -f 5`
            fasta_file=`echo $filename | cut -d'.' -f1`
            fasta_file=$fasta_file".fasta"
            break
        done

        makeblastdb -dbtype nucl -in ./temp/raw_prophages/$fasta_file

        for file in ./temp/groupings/$name/*
        do
            r_or_no_r=""
            filename=`echo $file | cut -d'/' -f 5`
            query_names=`echo $filename | cut -d'.' -f1`
            query_fasta=$query_names".fasta"
            query_file="./temp/raw_prophages/"$query_fasta
            cat $query_file > temp_query.fasta
            blastn -db ./temp/raw_prophages/$fasta_file -query ./temp_query.fasta -outfmt 6 | sort -nrk4 > ./temp_blast_res.txt
            r_or_no_r=`head -1 ./temp_blast_res.txt | awk '$9 > $10 {print "R"}'`
            head -1 ./temp_blast_res.txt | awk '$9 > $10 {print $1}' >> ./temp/easyfig_rc_log.txt

            feature_files=$feature_files" ./output/gbf/"$filename" "$r_or_no_r

            declare -i LYS_genes=`grep "/colour=32 160 64" $file | wc -l`
            printf $filename"\t"$LYS_genes"\n" >> LYS_length_log.txt

            declare -i REG_genes=`grep "/colour=33 96 192" $file | wc -l`
            printf $filename"\t"$REG_genes"\n" >> REG_length_log.txt

            declare -i EFF_genes=`grep "/colour=224 128 128" $file | wc -l`
            printf $filename"\t"$EFF_genes"\n" >> EFF_length_log.txt

            declare -i REC_genes=`grep "/colour=166 202 240" $file | wc -l`
            printf $filename"\t"$REC_genes"\n" >> REC_length_log.txt

            declare -i STR_genes=`grep "/colour=160 224 128" $file | wc -l`
            printf $filename"\t"$STR_genes"\n" >> STR_length_log.txt

            declare -i TRN_genes=`grep "/colour=96 32 128" $file | wc -l`
            printf $filename"\t"$TRN_genes"\n" >> TRN_length_log.txt

            declare -i MET_genes=`grep "/colour=192 160 192" $file | wc -l`
            printf $filename"\t"$MET_genes"\n" >> MET_length_log.txt

            declare -i OTH_genes=`grep "/colour=192 192 192" $file | wc -l`
            printf $filename"\t"$MET_genes"\n" >> OTH_length_log.txt

            declare -i IS_genes=`grep "/colour=255 255 0" $file | wc -l`
            printf $filename"\t"$IS_genes"\n" >> IS_length_log.txt

        done

        rm ./temp/raw_prophages/*.fasta.n*
        rm ./temp_query.fasta
        rm ./temp_blast_res.txt

    else
        for file in ./temp/groupings/$name/*
        do
            filename=`echo $file | cut -d'/' -f 5`

            declare -i LYS_genes=`grep "/colour=32 160 64" $file | wc -l`
            printf $filename"\t"$LYS_genes"\n" >> LYS_length_log.txt

            declare -i REG_genes=`grep "/colour=33 96 192" $file | wc -l`
            printf $filename"\t"$REG_genes"\n" >> REG_length_log.txt

            declare -i EFF_genes=`grep "/colour=224 128 128" $file | wc -l`
            printf $filename"\t"$EFF_genes"\n" >> EFF_length_log.txt

            declare -i REC_genes=`grep "/colour=166 202 240" $file | wc -l`
            printf $filename"\t"$REC_genes"\n" >> REC_length_log.txt

            declare -i REP_genes=`grep "/colour=192 160 192" $file | wc -l`
            printf $filename"\t"$REP_genes"\n" >> REP_length_log.txt

            declare -i STR_genes=`grep "/colour=160 224 128" $file | wc -l`
            printf $filename"\t"$STR_genes"\n" >> STR_length_log.txt

            declare -i TRN_genes=`grep "/colour=96 32 128" $file | wc -l`
            printf $filename"\t"$TRN_genes"\n" >> TRN_length_log.txt

            declare -i MET_genes=`grep "/colour=192 160 192" $file | wc -l`
            printf $filename"\t"$MET_genes"\n" >> MET_length_log.txt

            declare -i OTH_genes=`grep "/colour=192 192 192" $file | wc -l`
            printf $filename"\t"$MET_genes"\n" >> OTH_length_log.txt

            declare -i IS_genes=`grep "/colour=255 255 0" $file | wc -l`
            printf $filename"\t"$IS_genes"\n" >> IS_length_log.txt

        done

        feature_files="./temp/groupings/"$group_name"/*"

    fi
    python $easyfig_path -f tRNA [rect] -f IS [frame] -f CDS [arrow] -svg -o ./output/easyfig_res/$group_name.svg $feature_files
done

#############################################################################
# - All the logs previously made are moved to their repective folder, and an#
# extra log is made looking only at IS proportion in clusters which have on #
# average more or equal to 1 IS (this log doesn't included clusters with 0  #
# IS elements.                                                              #
# - An R script is called to graph all these logs and their trends, as well #
# as test their statistical significance.                                   #
#############################################################################

rm ./temp/length_logs/* &> /dev/null
mv ./*length_log.txt ./temp/length_logs/

for file in ./temp/length_logs/*length_log.txt
do
    filename=`echo $file | cut -d'/' -f 4 | cut -d'.' -f1`
    grep "cluster*" $file > ./temp/length_logs/$filename.summary.tab
done

cat ./temp/length_logs/IS_length_log.summary.tab | awk 'int($4)>0' > ./temp/length_logs/IS_length_log.summary2.tab

Rscript ./background/scripts/plot_script.R
