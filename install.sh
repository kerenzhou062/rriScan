#/bin/sh

base=`pwd`

echo "start to install rriScan..."

cd $base/thirdUtils/BamTools
make clean_api
make api
message=`make`

cd $base/thirdUtils/cdflib
make clean
message=`make`

cd $base/thirdUtils/RNAfoldLib
make clean
message=`make`

cd $base/bioUtils
make clean
message=`make`

cd $base
make clean
message=`make`

cd $base

echo "installation finished!"
