# Change Log

## 2017-Mar-25: Remove GUI-related code
Remove guimod_stub.f90, and comment out all the GUI-related code in puma.f90 with a description "XW(Mar/25/2017) to remove GUI:". The number of source code were shinked from 7 to 6.

## 2017-Mar-23: Make PUMA running

Successfully make PUMA running on 1CPU, with fake MPI and GUI interface (stub). 
I changed "nresources" line in puma.f90, which might be related to retrieving system time used in pumax.c.
They are marked with "XW".
