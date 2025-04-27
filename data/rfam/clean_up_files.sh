find . -name RF*.seed -print0 | xargs -0 rm -f
find . -name *.cov -print0 | xargs -0 rm -f
find . -name *.helixcov -print0 | xargs -0 rm -f
find . -name *.power -print0 | xargs -0 rm -f
find . -name *.surv -print0 | xargs -0 rm -f
find . -name *.sto -print0 | xargs -0 rm -f
find . -name *.pdf -print0 | xargs -0 rm -f
find . -name *.svg -print0 | xargs -0 rm -f
find . -name *.done -print0 | xargs -0 rm -f
find . -name '*_rnartist.svg' -print0 | xargs -0 rm -f
find . -name '*.cm' ! -name 'RF00005.cm' ! -name 'all.cm' -print0 | xargs -0 rm -f
