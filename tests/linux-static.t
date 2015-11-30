
Make sure we got static linkage if on Linux.

  $ ADMIXTURE=$TESTDIR/../admixture
  $ if [ `uname` != Linux ]; then exit 80; fi  
  $ readelf -d $ADMIXTURE
  
  There is no dynamic section in this file.

 