#!/bin/bash
sim (){
	java $JAVA_LIB -Xmx7g rccto3d.optimisation.MultiSearch $1.fa $1.pdb 1000 > out/rep.out
}
ref (){
	java $JAVA_LIB -Xmx7g rccto3d.optimisation.Refiner out/ 5 > out/repRef.out
}

JAVA_LIB=-Djava.library.path=.
DIREC=testEfficiency/use/

#crea archivos fasta y corige errores de formato en pdb
#javac extractSeqPDB.java
#mkdir -p testEfficiency/use/
#for X in testEfficiency/*; do
#	java extractSeqPDB $X A testEfficiency/use/
#done
#java extractSeqPDB testEfficiency/2p1g.pdb B testEfficiency/use/
#rm testEfficiency/use/testSetCompEfficiency*
#rm extractSeqPDB.class

make gccSimul
make
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
		res2=$(date +%s.%N)
		echo "Elapsed:    $(echo "$res2 - $res1"|bc )" > out/time.out
		mv out/* ${DIREC}final/$XX 
		;;
  	(*);;
	esac
done

