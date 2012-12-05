#!/bin/sh
#/bin/rm -Rf Output tempDir

tempDir=/data/glx-herschel2/data1/abeelen/sanepicInternal/test

/bin/rm -Rf Output tempDir $tempDir

mkdir $tempDir
ln -s  $tempDir tempDir


sanePre sanepic.ini 
sanePos sanepic.ini 
saneInv sanepic.ini

/bin/rm -f timing
for i in `seq 10`; do 
    echo $i;
    echo $i >> timing;
    /usr/bin/time -avpo timing sanePic sanepic.ini;
done
