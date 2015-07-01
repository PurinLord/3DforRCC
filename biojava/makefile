default:
	javac -d . aaClass/* PDBfromFASTA.java RotChain.java Trans.java MakeRClass.java ReformatAAPDB.java ArrangeAAPDB.java

trans:
	javac -d . Trans.java

pdbfrom:
	javac -d . PDBfromFASTA.java

aaClass:
	javac -d . MakeRClass.java

reformat:
	javac  -d . ReformatAAPDB.java

arrange:
	javac  -d . ArrangeAAPDB.java

rotamer:
	javac  -d . RotamerManager.java

simul:
	javac -d . SimulatedAnnealing3DProtFromRCC.java

test:
	javac TestRot.java

clearAll:
	rm -v *.class rccto3d/*.class rccto3d/rotamers/*.class

clear:
	rm -v rccto3d/*.class rccto3d/rotamers/*.class

clearThis:
	rm -v *.class
