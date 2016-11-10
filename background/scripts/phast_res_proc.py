#!/usr/bin/python
import os.path
import sys
from subprocess import call
import shutil

#############################################################################
# - Open system argument as readable file, followed by opening the output   # 
# file. Variables are reset prior to the loop. First line of the file is    #
# stored in variable line A.                                                #
#############################################################################

with open(sys.argv[1]) as file:
    region_out = open("./temp/phast_extrac_res/temp/regions_out_temp.txt","w")
    regions = ""
    region_final = ["",""]
    line_A = file.readline()

#############################################################################
# - The loop starts, and continues as long as lines are TRUE. Line B will   #
# read the second line. As the loop progresses line B transfers its data to #
# line A as more lines are read into line B. This means there will always be#
# one line difference between line A and line B variables. The loop ends    #
# once line A returns as empty (end of the input file).                     #
# - The variable contain start and end coordinates of prophage regions.     #
# - The variables are split into start and end coordinates for every line.  #
# If line B is found empty then its coordinates are replaced with the ending#
# coordinate of line A.                                                     #
#############################################################################

    while True:
        line_B = file.readline()

        if line_A == "": break
        
        data_A = line_A.split("-")
        data_B = line_B.split("-")
        region_A = int(data_A[0])
        region_B = int(data_A[1])

        if line_B == "":
            region_C = region_B
            region_D = region_B
        else:
            region_C = int(data_B[0])
            region_D = int(data_B[1])

        region_list = [region_A, region_B, region_C, region_D]

#############################################################################
# - At the end of the loop the region coordinate are added to the regions   #
# variable and written to the output file.                                  #
# - The comparison of variable is done through nested IF loops. The inner IF#
# loops check to see if the current coordinates are already within the      #
# final regions variable, if that is the case that means that this prophage #
# was merged in the prior loop, and a message is printed to the console.    #
# - The outter IF loop first checks whether there are 4000 bp between the   #
# end and start of the prophages, if yes the coordinates are not merged and #
# recorded normally. The second part of the IF loop checks for distances    #
# between prophages which is smaller or equal to 0 (overlapping prophages)  #
# if they are overlapping the coordinates are actually recorded normally    #
# (this may seem weird as we are merging prophages that are close together  #
# but the idea here is that if PHAST results have 3 overlapping prophages   #
# most likely these are neighboring prophages and PHAST is having issues    #
# calling boundaries, therefore, these should be kept as individuals. This  #
# was observed with the Sakai prophages with SP6,7,8). Finally if these two #
# are not the case this means that the prophages are within 4000 bps of each#
# other without overlapping and will be merged.                             #
#############################################################################

        if (int(region_list[2])-int(region_list[1])) > 4000:
            if str(region_list[1]) in regions:
                print "Prophages merged."
            else:
                region_final[0] = region_list[0]
                region_final[1] = region_list[1]

        elif (int(region_list[2])-int(region_list[1])) <= 0:
            if str(region_list[1]) in regions:
                print "Prophages merged."
            else:
                region_final[0] = region_list[0]
                region_final[1] = region_list[1]

        else:
            if str(region_list[1]) in regions:
                print "Prophages merged."
                region_final[1] = region_list[3]
            else:
                region_final[0] = region_list[0]
                region_final[1] = region_list[3]

        regions = str(region_final[0])+":"+str(region_final[1])+"\n"
        region_out.write(regions)    
        line_A = line_B

region_out.close()
