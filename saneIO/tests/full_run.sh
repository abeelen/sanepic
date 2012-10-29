#!/bin/sh
/bin/rm -Rf Output tempDir
mpirun.mpich2 -n 2 sanePre_PARA_FRAME sanepic.ini 
mpirun.mpich2 -n 2 sanePos_PARA_FRAME sanepic.ini 
mpirun.mpich2 -n 2 saneInv_PARA_FRAME sanepic.ini 
mpirun.mpich2 -n 24 sanePic_PARA_BOLO sanepic.ini 
