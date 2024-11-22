# Molpro interface

This is the version 2022.07 of the Molpro Interface. 
It converts an surface in internal coordinates to a analytic form and allows to evaluate it at a given geometry

# Basic Usage

1. Run 'make' in your terminal to build the version. Use the executable on your desired Input file to convert the 
   grid representation to an analytical representation.
2. Copy the files 'src/extern/*.F90' into your program code. In the file 'mod_poly.F90', you will find the routines
   to evaluate the surface correctly. There, you also find a better description
3. Keep in mind, that there are LAPACK and BLAS routines used in the code.

# Input File

First, please give the name of your potfile with

   pot='potfile.pot'

There are additional keywords to adapt the fit. Those are optional:

For coordinates, specific options can be set with the following command scheme:

coord=x,key1=l,key2=n,key3=k,...

possible keywords:
fitfct = (gauss/bspline/poly/morse/trigo)

# TODOS
 -> Add information of the geometry to be inputed
 -> Put more Information like number of grid points, or even ranges in the outfile
 -> Allow limitation of coupling terms in routine
 -> Make Keywords to change fit functions
