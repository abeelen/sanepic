#!/bin/sh
/bin/rm -Rf Output tempDir
sanePre sanepic.ini 
sanePos sanepic.ini 
saneInv sanepic.ini

/bin/rm -f timing
for i in `seq 3`; do 
    echo $i;
    echo $i >> timing;
    /usr/bin/time -avpo timing sanePic sanepic.ini;
done
