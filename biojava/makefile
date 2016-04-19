default:
	javac -d . aaClass/* PDBfromFASTA.java Trans.java SimulatedAnnealing.java StructureMannager.java MultiSearch.java Substitutor.java Refiner.java

trans:
	javac -d . Trans.java

pdbfrom:
	javac -d . PDBfromFASTA.java

optim:
	javac -d . SimulatedAnnealing.java StructureMannager.java MultiSearch.java Substitutor.java Heater.java Refiner.java

simul:
	javac -d . SimulatedAnnealing.java

struc:
	javac -d . StructureMannager.java

gccSimul:
	javah -jni -d ../SPalignNS/ rccto3d.optimisation.SimulatedAnnealing3DProtFromRCC
	g++ -I/usr/lib/jvm/java-7-openjdk-amd64/include/ -I/usr/lib/jvm/java-7-openjdk-amd64/include/linux/ -o libSPalig.so -shared -fpic ../SPalignNS/sim2sp.cc ../SPalignNS/sim2spImpl.cpp ../SPalignNS/sim2spImpl.h ../SPalignNS/protein2.cc ../SPalignNS/protein2.h ../SPalignNS/align1.cc ../SPalignNS/quatfit.cc ../SPalignNS/sp_misc.h ../SPalignNS/alignNS.cc ../SPalignNS/rccto3d_optimisation_SimulatedAnnealing3DProtFromRCC.h ../SPalignNS/sp_type.h ../SPalignNS/spalign.cc  

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
