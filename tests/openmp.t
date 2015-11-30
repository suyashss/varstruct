
Make sure OpenMP gives the same results.

  $ ADMIXTURE=$TESTDIR/../admixture
  $ DATA=$TESTDIR/../Data
  $ GOLD=$TESTDIR/gold
  $ $ADMIXTURE -j2 $DATA/hapmap3.bed 3 > /dev/null
  $ echo $?
  0
  $ diff -u $GOLD/hapmap3.3.Q  hapmap3.3.Q
  $ diff -u $GOLD/hapmap3.3.P  hapmap3.3.P

