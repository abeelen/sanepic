#!/bin/sh


/bin/rm -Rf Output 

export PM="mpirun.openmpi -n 72 -hostfile /home/abeelen/.openmpi.herschel.host -host glx-herschel,glx-herschel2" # mixed
export PM="mpirun.openmpi -n 4  -hostfile /home/abeelen/.openmpi.herschel.host -host glx-herschel"  # frame
export PM="mpirun.openmpi -n 24 -hostfile /home/abeelen/.openmpi.herschel.host -host glx-herschel" #bolo


$PM saneFrameOrder_PARA sanepic.ini
# $PM saneCheck_PARA sanepic.ini
# $PM saneFix_PARA sanepic.ini

/usr/bin/time -o Output/time_sanePre -v $PM sanePre_PARA sanepic.ini 
/usr/bin/time -o Output/time_sanePos -v $PM sanePos_PARA sanepic.ini 
/usr/bin/time -o Output/time_saneInv -v $PM saneInv_PARA sanepic.ini 
/usr/bin/time -o Output/time_sanePic -v $PM sanePic_PARA sanepic.ini 
