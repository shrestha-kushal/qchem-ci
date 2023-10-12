#!/bin/bash
ifort -O3 -c qchem_ci.f90 
ifort -O3 -c simplex.f
ifort -openmp -mkl=parallel -O3 simplex.o qchem_ci.o -lpthread -lm -o qchem_ci.exe
