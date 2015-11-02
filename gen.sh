#!/bin/bash
# first argument is type of dataset we're going to generate
E_BADARGS=65

if [ $# -ne 3 ]
then
    echo "Usage: `basename $0` {cmd} {type of dataset} {size}"
    exit $E_BADARGS
fi
SET=$1
SIZE=$3
TYPE=$2
FILE="${SET}-t${TYPE}-${SIZE}n.arff"
make

./"${SET}" -n $SIZE -t $TYPE
DATA="${SET}-t${TYPE}-${SIZE}n.dat"
if [[ ! -f $DATA ]]; then
  echo "data set ${DATA} was not found"
  exit 1
fi
tab2csv $DATA $FILE

cat <<EOF > "/tmp/out"
@RELATION ${set}-t${TYPE}-${SIZE}n

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
