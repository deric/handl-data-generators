#!/bin/csh
#
set echo
#
cp random_data.H /$HOME/include
#
g++ -c -g random_data.C >& compiler.out
if ( $status != 0 ) then
  echo "Errors compiling random_data.C."
  exit
endif
rm compiler.out
#
mv random_data.o ~/lib/$ARCH/random_data.o
#
echo "A new version of random_data has been created."
