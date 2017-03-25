# Change Log

## 2017-Mar-25: Add Post-Processing (pp)

Add directory "pp" for the post-processing excutable "burn7.x" and 2 associated namelist, one for 3 surface variables, another for 10 multi-level variables. To compile burn7.x you should apt install libnetcdf-cxx-legacy-dev first.

## 2017-Mar-25: Remove GUI-related code
Remove guimod_stub.f90, and comment out all the GUI-related code in puma.f90 with a description "XW(Mar/25/2017) to remove GUI:". Please note I did not delete any lines, just comment out. The number of source code were shinked from 7 to 6.

## 2017-Mar-23: Make PUMA running

Successfully make PUMA running on 1CPU, with fake MPI and GUI interface (stub). 
I changed "nresources" line in puma.f90, which might be related to retrieving system time used in pumax.c.
They are marked with "XW".
