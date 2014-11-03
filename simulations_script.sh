#!/bin/sh

# This file explains how to replicate the analyses performed in the 
# manuscript "Flexible methods for estimating genetic distances
# from nucleotide data" by Joly, Bryant and Lockhart.

# Go in folder src and compile the analysis software "anal" using
# the following command:

cd src
make install
cd ..

# You might need to update your gcc version to gcc 4.8 to properly compile
# as the program uses regex which are not fully implemented in previous versions
# of the compiler. A compiled version of the program is available for MAC
# intel64.

# For the remainning simulation, we expect that you have install fastsimcoal2.5
# (fsc25) in your path. You could also run it locally, but in such a case, you will
# have to add "./" before "fsc25" each time you want to run it.

################
#
# Figure 1 and 2
#

# Run fsc25 to perform the simulations used for figure 1:

fsc25 -i Fig1a.par -n 2000 --quiet --cores 2
fsc25 -i Fig1b.par -n 2000 --quiet --cores 2

# Then run "anal" to analyse the simulations output

./anal -i ./Fig1a/Fig1a_1.arb -c 1000
./anal -i ./Fig1b/Fig1b_1.arb -c 1000


###############################################
#
# Figure S1 and 3 - Rate variation among genes
#

fsc25 -t FigS1a.tpl -n 1 -f FigS1.def --quiet
fsc25 -t FigS1b.tpl -n 1 -f FigS1.def --quiet

# Now you need to concatenate all simulated file names into one document

cat ./FigS1a/*.arb > ./FigS1a/FigS1a_files.arb
cat ./FigS1b/*.arb > ./FigS1b/FigS1b_files.arb

# Then you need to analyse the results...

./anal -i ./FigS1a/FigS1a_files.arb -c 1000
./anal -i ./FigS1b/FigS1b_files.arb -c 1000


##################################################
#
# Figure S2 - Recombination (short length markers)
#

# Run fsc25 using the following command

fsc25 -i recomb.short.a.par -n 2000 --quiet --cores 2
fsc25 -i recomb.short.b.par -n 2000 --quiet --cores 2

# Then, run anal with the following command

./anal -i ./recomb.short.a/recomb.short.a_1.arb -c 1000
./anal -i ./recomb.short.b/recomb.short.b_1.arb -c 1000


##########################
#
# Figure 4 - Recombination
#

# Run fsc25 using the following command

fsc25 -i recomb.none.par -n 2000 --quiet --cores 2
fsc25 -i recomb.high.par -n 2000 --quiet --cores 2

# Then, run anal with the following command

./anal -i ./recomb.none/recomb.none_1.arb -c 10000
./anal -i ./recomb.high/recomb.high_1.arb -c 10000


###############################
#
# Figure 5 - Hybrid simulations
#

# Run fsc25 using the following command

fsc25 -i Fig5a.young.par -n 2000 --quiet --cores 1
fsc25 -i Fig5a.mid.par -n 2000 --quiet --cores 1
fsc25 -i Fig5a.old.par -n 2000 --quiet --cores 1

# Then run "anal" to analyse the simulations output

./anal -i ./Fig5a.young/Fig5a.young_1.arb -c 1000 -h
./anal -i ./Fig5a.mid/Fig5a.mid_1.arb -c 1000 -h
./anal -i ./Fig5a.old/Fig5a.old_1.arb -c 1000 -h

# Then rerun "anal" to analyse the simulations output
# assuming unequal contribution of the parents

./anal -i ./Fig5a.young/Fig5a.young_1.arb -c 1000 -h 1 -o Fig5b.young_results.txt
./anal -i ./Fig5a.mid/Fig5a.mid_1.arb -c 1000 -h 1 -o Fig5b.mid_results.txt
./anal -i ./Fig5a.old/Fig5a.old_1.arb -c 1000 -h 1 -o Fig5b.old_results.txt


########################
#
# Figure generation in R
#

# Open the script "Script_figures.R" in R and execute it to reproduce the
# figures of the manuscript. You can also execute the script using command line:

R CMD BATCH script_figures.R
