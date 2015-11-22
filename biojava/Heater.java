package rccto3d.optimisation;

import java.lang.Math;
import java.lang.Float;
import java.io.*;
import java.io.RandomAccessFile;

import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.AminoAcid;
import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.StructureTools;

import org.biojava.bio.structure.align.*;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.align.fatcat.calc.FatCatParameters ;
import org.biojava.bio.structure.align.fatcat.FatCatRigid;

import rccto3d.*;
import rccto3d.Trans;
import rccto3d.PDBfromFASTA;
import rccto3d.optimisation.SimulatedAnnealing3DProtFromRCC;

public class Heater{

public static void main(String args[]){
	
	Structure struc_ini = PDBfromFASTA.readPDB(args[0]);
	Structure struc_fin = PDBfromFASTA.readPDB(args[0]);
	float minRMSd = Float.parseFloat(args[1]);
	float variation = minRMSd/100;
	double currentRMSd = 0.0;
	int searchCount = 0;
	
	SimulatedAnnealing3DProtFromRCC simA = new SimulatedAnnealing3DProtFromRCC("", "");
	while(currentRMSd < minRMSd){
		struc_fin = simA.alterConformationAll(struc_fin, variation, variation);
		
		currentRMSd = simA.calcRMSD2(struc_ini, struc_fin);
		System.out.println("RM2> " + currentRMSd);
		//currentRMSd = simA.calcRMSD(struc_ini, struc_fin);
		//System.out.println("RMS: " + currentRMSd);

		searchCount ++;
	}

	//foma MUY idiota de guardar en directorio MUY
	String[] dirs = args[0].split("/");
	int last = dirs.length;
	String name = "";
	String dir = "";
	for (String i : dirs){
		if(i.substring(i.length() - 4).equals(".pdb")){
			name = i.substring(0, i.length() - 4);
		}else{
			dir = dir + i + "/";
		}

	}
	//System.out.println(dir + " " + name);
	System.out.println("RMSd obtenido: " + currentRMSd);
	System.out.println(("variación promedio: " + variation/2 * searchCount));
	//struc_fin.setName("RMSd obtenido: " + currentRMSd +"\tvariación promedio: " + variation/2 * searchCount);
	String fileName = dir + name + "_" + args[1] + ".pdb";
	PDBfromFASTA.writePDB(fileName, struc_fin);
}
}
