
These are some simple tests using the example dataset.

  $ VARSTRUCT=$TESTDIR/../varstruct
  $ DATA=$TESTDIR/../Data
  $ $VARSTRUCT -d $DATA/exampledata.stru -k 1 -n 60 -m 100 -p 2 -o testout > /dev/null
  $ echo $?
  0
  $ $VARSTRUCT -d $DATA/exampledata.stru -k 2 -n 60 -m 100 -p 2 -o testout > /dev/null
  $ echo $?
  0
  $ $VARSTRUCT -d $DATA/exampledata.stru -k 3 -n 60 -m 100 -p 2 -o testout > /dev/null
  $ echo $?
  0
  $ $VARSTRUCT -d $DATA/exampledata.stru -k 3 -n 60 -m 100 -p 2 -o testout > /dev/null
  $ echo $?
  0

