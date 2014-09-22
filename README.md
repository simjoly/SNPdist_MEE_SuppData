SNPdist_MEE_SuppData
====================

This folder contains all the scripts and code necessary to replicate the analyses performed in the manuscript "Flexible methods for estimating genetic distances from nucleotide data" by Joly, Bryant and Lockhart published in MEE.

The present file explains how to replicate the analyses.


Instructions
------------

All the commands needed to replicate the simulations and analyses performed in the manuscript can be found in the shell script "simulations_script.sh". Therefore, you could just type 

$> sh simulations_script.sh

to perform all the analyses. However, the script is well commented and you can explore it to run the analysis one step at a time.


Requirements
------------

h3. fsc25

We expect that you have installed fastsimcoal version 2.5 (fsc25; version 2.5.0.2 was used in the paper) and that you have placed it in your $PATH. You could also run it locally, but in such a case, you will have to add "./" before "fsc25" each time you want to run it (and modify the shell script accordingly).

h3. gcc

You also need to have a version of the gcc compiler installed on your computer to compile the program anal that analyses the output of fastsimcoal. You might need to update your gcc version to gcc 4.8 to properly compile "anal" as the program uses regex which is not fully implemented in previous versions of the compiler. A compiled version of the program is available for MAC intel64.

h3. R

If you want to reproduce the figures used in the manuscript, you will also need to have R installed on you computer, with the librarie 'plyr', 'ggplot', and 'gridExtra'.
