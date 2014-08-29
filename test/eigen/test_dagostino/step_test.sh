#!/bin/sh

usage(){
cat << EOF

ERROR 

usage : ./step_test.sh np start input

np : number of processors used
start : number of node of the matrix to start with
input : file containing the entry of test_realnet

EOF
}

if test $# -ne 3; then
        usage
        exit 0
fi
file=$3
start=$2
np=$1

#newfile="$file"new
#sed -i 1d enea.inp
#sed -i '1i'$newfile $n enea.inp

#cp $file $newfile
n=`expr $start`
sed -i 1d $file
sed -i '1i'$n $file
echo "on rentre dans la boucle"
while [ $n -le 7799880 ] ;do
  	sed -i 1d $file
	sed -i '1i'$n $file
	echo "fichier modifié"
	echo $n
	echo
	mpirun -np $np ./test_realnet < $file
	n=$(($n+10000))
done
