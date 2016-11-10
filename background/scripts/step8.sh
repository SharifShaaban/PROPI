#############################################################################
# - The average proportion of prophage to the strain's chromosome is        #
# calculated by adding the lengths of all the prophages of one strain, and  #
# getting the percentage that represents of the strain's chromosome.        #
# - Tricks like "rev" are used to quickly isolate the prophage length       #
# without having to go through the file.                                    #
#############################################################################

rm ./output/p_to_c_proportions.txt &> /dev/null

for genome in ./input/sequences/*.fasta
do
    g_length=`tail -n +2 "$genome" | wc -c`
    g_name=`echo $genome | cut -d'/' -f 4 | cut -d'_' -f 1`
    declare -i int_glength=$g_length
    declare -i int_plength=0
    for prophage in ./output/gbf/$g_name*.gbf
    do
        p_length=`echo $prophage | cut -d'/' -f 4 | cut -d'.' -f 1 | rev | cut -d'_' -f1 | rev`
        declare -i int_splength=$p_length
        int_plength=$int_plength+$int_splength
    done
    printf "\n"$g_name"\n" >> ./output/p_to_c_proportions.txt
    proportion=`echo "scale=2; $int_plength * 100 / $int_glength" | bc -l`
    printf $int_glength" "$int_plength" "$proportion"\n" >> ./output/p_to_c_proportions.txt
done
