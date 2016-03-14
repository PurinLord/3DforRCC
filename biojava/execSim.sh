#!/bin/bash
sim (){
	java $JAVA_LIB rccto3d.optimisation.MultiSearch $1.fa $1.pdb 1000 > out/rep.out
}
ref (){
	java $JAVA_LIB rccto3d.optimisation.Refiner out/ 5 > out/repRef.out
}

JAVA_LIB=-Djava.library.path=.
DIREC=testEfficiency/use/
echo ${DIREC}
#for X in ${DIREC}*.pdb; do
mkdir -p ${DIREC}final/
for X in `ls ${DIREC}`; do
	case $X in
  	(*.pdb)
		res1=$(date +%s.%N)
		XX=`echo $X | sed 's/\(.*\)\.pdb/\1/'`
		echo $XX
		mkdir -p ${DIREC}final/$XX 
		sim $DIREC$XX
		ref $DIREC$XX
		#(java rccto3d.optimisation.MultiSearch $XX.fa $XX.pdb 1000 > out/rep.out);(java rccto3d.optimisation.Refiner out/ 10 > repRef.out); 
		res2=$(date +%s.%N)
		echo "Elapsed:    $(echo "$res2 - $res1"|bc )" > out/time.out
		mv out/* ${DIREC}final/$XX 
		;;
  	(*);;
	esac
done

