#############################################################################
# - Run GetHomologues at different coverage and similarity percentages      #
# accross all the isolated prophage sequences.                              #
#############################################################################

rm -r ./temp/*_homologues/

for per in `cat ./input/percent_list.txt`
do
    per_com="_C"$per"_S"$per
    per_com_dir=$per_com"_"
    cd ./temp/

    /groups2/gally_grp/fsa_project/Nadya/software/get_homologues-x86_64-20140930/get_homologues.pl -c -n 8 -e -C $per -S $per -d ./temp_gbf/
    /groups2/gally_grp/fsa_project/Nadya/software/get_homologues-x86_64-20140930/get_homologues.pl -M -n 8 -e -t 0 -C $per -S $per -d ./temp_gbf/
    /groups2/gally_grp/fsa_project/Nadya/software/get_homologues-x86_64-20140930/get_homologues.pl -G -n 8 -e -t 0 -C $per -S $per -d ./temp_gbf/
# Repeat of gethomologues step as it tends to not finish on first run when too many samples are used, but finishes second run correctly
    /groups2/gally_grp/fsa_project/Nadya/software/get_homologues-x86_64-20140930/get_homologues.pl -M -n 8 -e -t 0 -C $per -S $per -d ./temp_gbf/
    /groups2/gally_grp/fsa_project/Nadya/software/get_homologues-x86_64-20140930/get_homologues.pl -G -n 8 -e -t 0 -C $per -S $per -d ./temp_gbf/
    cd ../

    bdbh_file=`ls -d ./temp/temp_gbf_homologues/*algBDBH*$per_com_dir/`
    omcl_file=`ls -d ./temp/temp_gbf_homologues/*algOMCL*$per_com_dir/`
    cog_file=`ls -d ./temp/temp_gbf_homologues/*algCOG*$per_com_dir/`
    /groups2/gally_grp/fsa_project/Nadya/software/get_homologues-x86_64-20140930/compare_clusters.pl -o ./core$per_com/ -m -T -n -d $bdbh_file, $cog_file, $omcl_file
    /groups2/gally_grp/fsa_project/Nadya/software/get_homologues-x86_64-20140930/compare_clusters.pl -o ./acc$per_com/ -m -T -n -d  $cog_file, $omcl_file
done

#############################################################################
# - Move GetHomolgues folder away to main results folder, and extract matrix#
# and tree moving them to dedicated folders.                                #
#############################################################################

mkdir ./output/gethomologues_res &> /dev/null
mkdir ./output/pangenome_trees &> /dev/null
mkdir ./temp/pangenome_matrix_bin &> /dev/null
rm -r ./output/gethomologues_res/*
mv ./core*/ ./output/gethomologues_res/
mv ./acc*/ ./output/gethomologues_res/

for folder in ./output/gethomologues_res/acc_C*
do
    foldername=`echo $folder | cut -d'/' -f 4`
    find $folder/ -name "*.tab" -exec cp {} ./temp/pangenome_matrix_bin/$foldername.pangenome_matrix.tab \;
done

#need to automate cp
cp ./output/gethomologues_res/acc_C75_S75/pangenome_matrix_t0.phylip.ph ./output/pangenome_trees/C75.ph

# removing _gbf... out of tree labels
for file in ./output/pangenome_trees/*.ph
do
    sed -i "s/_gbf...//g" $file
done
