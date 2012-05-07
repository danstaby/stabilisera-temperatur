#! /bin/bash

ERROR="Too few arguments : no file name specified"
[[ $# -eq 0 ]] && echo $ERROR && exit # no args? ... print error and exit

# check that the file exists
if [ -f $1.tex ]
then
# if it exists then latex it twice, dvips, then ps2pdf, then remove all the unneeded files
latex $1.tex
bibtex $1.aux
latex $1.tex
dvips $1.dvi -o $1.ps
ps2pdf $1.ps

else
# otherwise give this output line with a list of available tex files
echo no such file! Choose one of these:
fi
