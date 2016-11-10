#!/usr/bin/python
import os.path

#############################################################################
# - Keywords are distributed into their respective functional list.         #
# - A list of RGB colours is made.                                          #
# - Outfiles are called in appending format.                                #
#############################################################################

regulation_list=["anti-repressor","anti-termination","antirepressor","antitermination","antiterminator","anti-terminator","cold-shock","cold shock","heat-shock","heat shock","regulation","regulator","regulatory","helicase","antibiotic resistance","repressor","zinc","sensor","dipeptidase","deacetylase","5-dehydrogenase","glucosamine kinase","glucosamine-kinase","dna-binding","dna binding","methylase","sulfurtransferase","acetyltransferase","control","ATP-binding","ATP binding","Cro","Ren protein","CII","inhibitor","activator","derepression","protein Sxy","sensing","sensor","Tir chaperone","Tir-cytoskeleton","Tir cytoskeleton","Tir protein","EspD"]
eff_vf_list=["effector","T3S","virulence","viral","Type III","Intimin","espA","multidrug","tape measure","toxin","Eae protein","EscD","EspF","enterohemolysin", "Lom protein","NleB","NleA8","secreted","secretion","SepL protein"]
rec_rep_list=["excisionase","integrase","transposase","resolvase","recombination","recombinase","recombinatory","replication","replicative","Insertion element","invertase","ninB","ninG","ninH","restriction","Kil protein","host specificity","IS","DnaB","DicF","DNA primase","DNA polymerase","transposition","protein RecT","RecT protein","TnpA","Cl protein","cspa","cytochrome"]
met_trans_list=["transport","metabolism","hydrolase","synthesis","glycosidase","synthase","synthetase","ribosomal","thiolase","lyase","pump","aminotransferase","amino-transferase","protease","phosphate reductase","uptake","periplasmic","sortase","permease","Bifunctional Protein PutA","Bifunctional protein TrpGD","kinase"]
structure_list=["scaffold","capsid","terminase","structural","tail","assembly","head","structure","fimbrial","flagellar","baseplate","filament","lipid","lipoprotein"]
lysis_list=["holin","lysozyme","DNA packaging","DNA-packaging","portal","DNA injection","DNA transfer","lysis","release","racemase","endopeptidase","carboxypeptidase","host killing","endolysin","endo-lysin","lisogenization","lytic"]
trna_list=["tRNA"]
stx_list=["shiga"]
colour_list=["33 96 192","224 128 128","166 202 240","192 160 192","160 224 128","32 160 64","224 32 0","96 32 128","192 192 192","224 96 0"]

out_reg = open("./temp/annot_groups/reg.tmp","a")
out_eff = open("./temp/annot_groups/eff.tmp","a")
out_rec = open("./temp/annot_groups/rec.tmp","a")
out_str = open("./temp/annot_groups/str.tmp","a")
out_lys = open("./temp/annot_groups/lys.tmp","a")
out_trn = open("./temp/annot_groups/trn.tmp","a")
out_stx = open("./temp/annot_groups/stx.tmp","a")
out_met = open("./temp/annot_groups/met.tmp","a")
out_oth = open("./temp/annot_groups/oth.tmp","a")

#############################################################################
# - For each file in the temp_gbf folder open, read the file line by line if#
# line contains the pattern "/product", determine in which function group it#
# belongs and create a colour flag line. In other words the gbf files are   #
# rewritten line by line with a colour flag line added after product lines. #
# - Products are also added to a temp file relevant to the functional group #
# it belongs to.                                                            #
#############################################################################

for filename in os.listdir("./temp/temp_gbf/"):

    with open("./temp/temp_gbf/"+filename) as file:
        line_A = file.readline()
        out_gbk = open(("./output/gbf/"+filename),"w")

        while True:
            if line_A == "//\n":
                out_gbk.write(line_A)
                break

            if "/product" in line_A:
                out_gbk.write(line_A)

                if any(stx_element.lower() in line_A.lower() for stx_element in stx_list):
                    out_gbk.write("                     /colour="+colour_list[6]+"\n")
                    out_stx.write(line_A)
                elif any(reg_element.lower() in line_A.lower() for reg_element in regulation_list):
                    out_gbk.write("                     /colour="+colour_list[0]+"\n")
                    out_reg.write(line_A)
                elif any(eff_element.lower() in line_A.lower() for eff_element in eff_vf_list):
                    out_gbk.write("                     /colour="+colour_list[1]+"\n")
                    out_eff.write(line_A)
                elif any(rec_element.lower() in line_A.lower() for rec_element in rec_rep_list):
                    out_gbk.write("                     /colour="+colour_list[2]+"\n")
                    out_rec.write(line_A)
                elif any(structure_element.lower() in line_A.lower() for structure_element in structure_list):
                    out_gbk.write("                     /colour="+colour_list[4]+"\n")
                    out_str.write(line_A)
                elif any(lys_element.lower() in line_A.lower() for lys_element in lysis_list):
                    out_gbk.write("                     /colour="+colour_list[5]+"\n")
                    out_lys.write(line_A)
                elif any(trn_element.lower() in line_A.lower() for trn_element in trna_list):
                    out_gbk.write("                     /colour="+colour_list[7]+"\n")
                    out_trn.write(line_A)
                elif any(met_element.lower() in line_A.lower() for met_element in met_trans_list):
                    out_gbk.write("                     /colour="+colour_list[3]+"\n")
                    out_met.write(line_A)
                else:
                    out_gbk.write("                     /colour="+colour_list[8]+"\n")
                    out_oth.write(line_A)

            else:
                out_gbk.write(line_A)
            
            line_A = file.readline()

exit
