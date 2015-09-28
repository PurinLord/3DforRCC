package rccto3d.optimisation;

import java.lang.Math;
import java.lang.Integer;
import java.io.*;

import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.AminoAcid;

import rccto3d.*;
import rccto3d.Trans;
import rccto3d.PDBfromFASTA;
import rccto3d.optimisation.SimulatedAnnealing3DProtFromRCC;

public class Heater{

public static void main(String args[]){
	
	Structure struc_ini = PDBfromFASTA.readPDB(args[0]);
	int heat = Integer.parseInt(args[1]);
	double variation = 10.0;
	
	for(; heat > 1; heat--){
		SimulatedAnnealing3DProtFromRCC.alterConformationAll(struc_ini, variation, variation);
	}

	String[] dirs = args[0].split("/");
	int last = dirs.length;
	String name = "";
	String dir = "";
	for (String i : dirs){
		if(i.substring(i.length() - 4).equals(".pdb")){
			name = i;
		}else{
			dir = dir + i + "/";
		}

	}
	//System.out.println(dir + " " + name);
	PDBfromFASTA.writePDB(dir + name + args[1] + "_" + variation , struc_ini);
}

}
