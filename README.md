# PROPI version 0.8.0

**PR**ophage **O**verview **PI**peline

Extracts prophage regions from sequences and clusters them based on gene content.

Required Dependencies:
-----------------------

* Python 2.6.6 (required modules: os.path, sys, subprocess, shutil)
* Python 2.7.9 for step7.sh (same required modules) (path can be changed in step7.sh)
* EMBOSS 6.5.7.0
* BLAST 2.2.28+
* R 3.0.0 (required libraries: magrittr, readr, ggplot2, cowplot)
* R 3.2.2 for step7.sh (same required libraries) (path can be changed in step7.sh)

* Perl 5.18.1 called in step4.sh line 7
* Prokka 1.5.2 called in step4.sh line 21 (Prokka requires its own dependencies here called in step4.sh line 8)
* Get_Homologues 1.0 called in step5.sh lines 15-20 and 26-27 (Get_Homologues requires its own dependencies)
* Easyfig_CL_2.1.py called in step7.sh line 3 (Easy_fig requires its own Python module dependencies)

**THE PATH FOR THESE TOOLS NEED TO BE CHANGED WITHIN THE SCRIPTS THEY ARE CALLED IN!!**

Required Inputs:
-----------------

* Sequences of interest (required to be above 100kbp) in the ./input/sequences/ folder
* A list of the sequences to be fully processed in the ./input/sequence_list.txt file (example of file format included)
* A list of BLAST identity and coverage percentage at which you want Get_Homologues to work (examle of file format included)

Running Pipeline:
------------------

* All the following commands should be run from the tool directory containing the background, temp, output, and input folders.
* Once every command runs it should create the required path, and shouldn't require any human interaction if given all correct permissions, dependencies, and inputs.
* The following commands are sequential and should be run once the previous command is fully completed.

Commands:
----------

``$ bash ./background/scripts/step1.sh``

``$ bash ./background/scripts/step2.sh``

``$ bash ./background/scripts/step3.sh``

``$ bash ./background/scripts/step4.sh``

``$ bash ./background/scripts/step5.sh``

``$ bash ./background/scripts/step6.sh``

``$ bash ./background/scripts/step7.sh``

``$ bash ./background/scripts/step8.sh``

Outputs (found in the ./output/ folder):
-----------------------------------------

* gethomologues_res folder: the results files from Get_Homologues looking at the core (if any) and accessory genomes.
* pangenome_trees folder: the outputted trees of Get_Homologues.
* gbf folder: all the outputted gbf files for prophages including functional colour flags for genes.
* prophage_code.tab: a Tab delimited file clustering the prophages by Eucledian distances.
* easyfig_res: easyfig alignments for all the t4.5 clusters detailed in prophage_code.tab.  
* easyfig_rc_log.txt: a log of all the prophages which had to be reverse complemented for the easyfig figures.
* gene_freq_plot_ggplot2.svg: an SVG plot outputted by R looking at prophage gene content based on gene function.

Comments:
----------

This pipeline can easily be modified to target specific areas of interest, comments are included in the scripts for guidance. If any questions please contact Sharif.Shaaban@roslin.ed.ac.uk.

This pipeline is set up to use the PHAST online tool. However, a newer version is now available PHASTER: a better, faster version of the PHAST phage search tool. Nucleic Acids Res., 2016 May 3) but has not yet been implemented in the pipeline.

Step4.sh contains an area that is fully commented out. This area was to obtain gene IDs for a DAVID analysis but requires annotations from a different source as PROKKA does not provide these type of IDs.

References:
------------

1- CAMACHO, C., COULOURIS, G., AVAGYAN, V., MA, N., PAPADOPOULOS, J., BEALER, K. & MADDEN, T. L. 2009. BLAST+: architecture and applications. BMC Bioinformatics, 10, 421.

2- CONTRERAS-MOREIRA, B. & VINUEASA, P. 2013. GET_HOMOLOGUES, a Versatile Software Package for Scalable and Robust Microbial Pangenome Analysis. Applied and Environmental Microbiology, 79.

3- RICE, P., LONGDEN, I. & BLEABSY, A. 2000. EMBOSS: The European Molecular Biology Open Software Suite (2000) Trends in Genetics, 16.

4- SEEMANN, T. 2014. Prokka: rapid prokaryotic genome annotation. Bioinformatics, 30, 2068-9.

5- SIGUIER, P., PEROCHON, J., LESTRADE, L., MAHILLON, J. & CHANDLER, M. 2006. ISfinder: the reference centre for bacterial insertion sequences. Nucleic Acids Res, 34, D32-6.

6- SULLIVAN, M. J., PETTY, N. K. & BEATSON, S. A. 2011. Easyfig: a genome comparison visualizer. Bioinformatics, 27, 1009-10.

7- WICKHAM, H. 2009. ggplot2: Elegant Graphics for Data Analysis., Springer-Verlag New York.

8- ZHOU, Y., LIANG, Y., LYNCH, K. H., DENNIS, J. J. & WISHART, D. S. 2011. PHAST: a fast phage search tool. Nucleic Acids Res, 39, W347-52.
