


dirs="nr iniparser saneIO saneLib saneFrameOrder sanePos sanePre sanePic saneInv sanePS saneFix saneCheck saneSplit saneMerge"
trash='*~ *.o *.la *.lo autoscan.log aclocal.m4 config.guess config.h config.sub install-sh ltmain.sh configure missing config.status Makefile.in Makefile config.log config.h.in libtool stamp-h1 depcomp autom4te.cache m4/libtool.m4 m4/lt~obsolete.m4 m4/ltoptions.m4 m4/ltsugar.m4 m4/ltversion.m4 .deps .libs'

rm -Rf ${trash};


for dir in ${dirs}; do
    echo ${dir};
    (cd ${dir}; rm -Rf ${trash});
    (cd ${dir}/src; rm -Rf ${trash});
    rm -f ${dir}/Release*
    rm -f ${dir}/Debug*
    rm -f ${dir}/src/${dir};
done

# autoreconf -i
