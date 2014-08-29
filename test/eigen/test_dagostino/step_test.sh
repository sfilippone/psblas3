#!/bin/sh

usage(){
cat << EOF

ERROR 

usage : ./step_test.sh np start totNzero input

np : number of processors used
start : number of matrix non zeros to start with
totNzero : max non zero entry's of the input matrix
input : file containing the input of test_realnet

EOF
}
if [ $1 = --help ]; then
	usage
	exit 0
fi

if test $# -ne 4; then
        usage
        exit 0
fi
file=$4
totNzero=$3
start=$2
np=$1

n=`expr $start`
sed -i '1i'$n $file
echo "on rentre dans la boucle"
while [ $n -le `expr $totNzero` ] ;do
  	sed -i 1d $file
	sed -i '1i'$n $file
	echo "fichier modifié"
	echo $n
	echo
	mpirun -np $np ./test_realnet < $file
	n=$(($n+10000))
done
sed -i 1d $file
