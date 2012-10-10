========================================================================

OBTNWP 0.4 written in Dept. of Atmos Sci, NJU, in 2004. 

========================================================================

OBTNWP 0.5 modify the I/O mode to 'stream', for compiler 
     and OS independant.

/////////////////////////////////////////////////////////////////////////////
Other notes:

Usage: 
  1. gfortran OBTNWP-0.5.utf8.f90 -o OBTNWP-0.5.utf8
  2. ./OBTNWP-0.5.utf8
  3. grads -blc "date.gs"
  4. gxeps -c -i Height.gmf -o Height.eps


To be added:
  1. Namelist replace .ini control file
  2. Makefile for compiling
  3. Adding NetCDF I/O
  4. Evolution to Baroclinic
 
/////////////////////////////////////////////////////////////////////////////
