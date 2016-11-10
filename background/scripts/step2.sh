#############################################################################
# - For the strains listed in the ./sequence_list.txt file send sequence    #
# in PHAST format to the PHAST server, then from the obtained file extract  #
# the URL for the summary of the results                                    #
# - Until the file is not of size 0 sleep for 5 minutes, then try accessing #
# the result URL again (until the result URL is populated by PHAST it will  #
# result in no file, thus the loop will check for results every 5 minutes.  #
# - Awk and grep are used to extract the coordinates of the PHAST results.  #
#############################################################################

phast_path=./temp/phast_extrac_res
mkdir $phast_path/html_out &> /dev/null
mkdir $phast_path/summary_out &> /dev/null
mkdir $phast_path/image_out &> /dev/null
mkdir $phast_path/temp &> /dev/null

for strain_name in `cat ./input/sequence_list.txt`
do
    wget --post-file $phast_path/seq_phast_fmt/$strain_name* http://phast.wishartlab.com/cgi-bin/phage_command_line.cgi -O $phast_path/html_out/$strain_name.phast_out.html
    url_res=`cat $phast_path/html_out/$strain_name.phast_out.html | sed 's/br>/\n/g' | grep http | sed 's/<//g' | grep summary`

    until [ -s $phast_path/summary_out/$strain_name.summary.txt ];
        do
            sleep 300
            wget $url_res -O $phast_path/summary_out/$strain_name.summary.txt
        done

    awk '{print $5}' $phast_path/summary_out/$strain_name.summary.txt  | grep "\-" > $phast_path/temp/$strain_name.temp
done

#############################################################################
# - Extract the URL for the PHAST image result and record it                #
#############################################################################

for file in $phast_path/html_out/*
do
    strain_name=`echo $file | cut -d"/" -f4 | cut -d"." -f2`
    image_url=`cat $file | sed 's/br>/\n/g' | grep http | sed 's/<//g' | grep image`
    wget $image_url -O $phast_path/image_out/$strain_name.png
done
