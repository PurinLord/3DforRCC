javac extractSeqPDB.java
mkdir testEfficiency/use/
for X in testEfficiency/*; do
	java extractSeqPDB $X A testEfficiency/use/
done

