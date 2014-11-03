#
# Linux/UNIX Makefile for Analysis version 1.00
#

# Objects
OBJS = information.o nexusdata.o organismsdata.o GenFunctions.o anal.o

# Flags
 CFLAGS = -O3 -std=c++11
# CFLAGS =

#If using gcc
#CC = gcc -fopenmp $(CFLAGS)

#If using cc
#CC = cc $(CFLAGS)

#If using g++
CC = g++ $(CFLAGS)

#Folder containning binairies of the program
EXEDIR = ../

a1b2c3d4:
	@echo " "
	@echo " "
	@echo "  To use the Analysis version 1.00 Makefile, type"
	@echo " "
	@echo "     make install       to compile Analysis vers. 1.00 and then install the "
	@echo "                          executables in $(EXEDIR), then remove the object"
	@echo "                          files to save space"
	@echo "     make all           to compile the whole package and leaving the executables"
	@echo "                          in the current directoty and leaving the object files"
	@echo "     make put           to move the executables into $(EXEDIR)"
	@echo "     make clean         to remove all object files and executables from the"
	@echo "                          current directory"
	@echo "     make pofad         to compile pofad and leaving the executables in the"
	@echo "                          current directoty and leaving the object files"
	@echo " "

introduce:
	@echo ""
	@echo "Building anal version 1.00"
	@echo ""

all:        introduce anal
	@echo ""
	@echo "Finished compiling"
	@echo ""

put:
	@echo "Moving anal version 1.00 binaries in $(EXEDIR)"
	@cp anal $(EXEDIR)

install:        all put clean
	@echo "Done compiling anal"

clean:
	@echo "Removing object files to save space"
	@rm -f *.o
	@rm -f anal


anal: $(OBJS)
	$(CC) $(OBJS) -o anal -lm

anal.o: anal.h 
	$(CC) -c anal.cpp

organismsdata.o: organismsdata.h
	$(CC) -c organismsdata.cpp 

nexusdata.o: nexusdata.h
	$(CC) -c nexusdata.cpp

information.o: information.h
	$(CC) -c information.cpp

GenFunctions.o: GenFunctions.h
	$(CC) -c GenFunctions.cpp

