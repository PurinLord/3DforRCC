default:
	javac -d . aaClass/* PDBfromFASTA.java RotChain.java Trans.java MakeRClass.java ReformatAAPDB.java ArrangeAAPDB.java

trans:
	javac -d . Trans.java

pdbfrom:
	javac -d . PDBfromFASTA.java

optim:
	javac -d . SimulatedAnnealing3DProtFromRCC.java MultiSearch.java Substitutor.java Heater.java Refiner.java

simul:
	javac -d . SimulatedAnnealing3DProtFromRCC.java

multi:
	javac -d . MultiSearch.java

ref:
	javac -d . Refiner.java

heat:
	javac -d . Heater.java

subst:
	javac -d . Substitutor.java

manager:
	javac -d . RCCManager.java

mrclass:
	javac -d . MakeRClass.java

reformat:
	javac  -d . ReformatAAPDB.java

arrange:
	javac  -d . ArrangeAAPDB.java

rotamer:
	javac  -d . RotamerManager.java

test:
	javac TestRot.java

clearAll:
	rm -v *.class rccto3d/*.class rccto3d/rotamers/*.class rccto3d/optimisation/*.class

clear:
	rm -v rccto3d/*.class rccto3d/rotamers/*.class

clearThis:
	rm -v *.class
