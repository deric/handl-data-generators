#!/bin/bash
# first argument is type of dataset we're going to generate
E_BADARGS=65

if [ $# -ne 2 ]
then
    echo "Usage: `basename $0` {type of dataset} {size}"
    exit $E_BADARGS
fi
SIZE=$2
TYPE=$1
FILE="cure-t${TYPE}-${SIZE}n-2D.arff"
make cure
./cure -n $SIZE -t $TYPE
DATA="cure-t${TYPE}-${SIZE}n-2D.dat"
if [[ ! -f $DATA ]]; then
  echo "data set ${DATA} was not found"
  exit 1
fi
tab2csv $DATA $FILE

cat <<EOF > "/tmp/out"
@RELATION cure-t${TYPE}-${SIZE}n-2D

@ATTRIBUTE x REAL
@ATTRIBUTE y REAL
@ATTRIBUTE class {0,1,2}

@DATA
EOF
cat "$FILE" >> /tmp/out && mv /tmp/out $FILE
if [[ -d ~/dev/clustering-benchmark ]]; then
  cp $FILE ~/dev/clustering-benchmark/src/main/resources/datasets/artificial/$FILE
  cd ~/dev/clustering-benchmark
  mvn install -DskipTests
fi
