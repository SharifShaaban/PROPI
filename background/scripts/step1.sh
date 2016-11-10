#############################################################################
# - For files in the sequence replace the first line "> ..." by "fasta_seq="#
# which is the required PHAST sequence format for bulk URLAPI submissions.  #
# - Use cut to obtain filename, then create new sequence in the PHAST       #
# sequence folder.                                                          #
#############################################################################

mkdir ./temp &> /dev/null
mkdir ./input &> /dev/null
mkdir ./input/sequences &> /dev/null
mkdir ./output &> /dev/null

phast_path=./temp/phast_extrac_res
mkdir $phast_path &> /dev/null

for strain in ./input/sequences/*
do
    filename=`echo $strain | cut -d'/' -f 4`
    sed '1 s/^.*$/fasta_seq=/g' $strain > $phast_path/seq_phast_fmt/$filename
done
