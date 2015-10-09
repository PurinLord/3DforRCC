package rccto3d.optimisation;

import java.lang.Math;
import java.lang.Integer;
import java.io.*;

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
	double minRMSd = Integer.parseInt(args[1]);
	double variation = 1.0;
	double currentRMSd = 0.0;
	
	while(currentRMSd < minRMSd){
		SimulatedAnnealing3DProtFromRCC.alterConformationAll(struc_fin, variation, variation);
		
		try {
		
			// To run FATCAT in the flexible variant say
			// FatCatFlexible.algorithmName below
			StructureAlignment algorithm  = StructureAlignmentFactory.getAlgorithm(FatCatRigid.algorithmName);

			Atom[] ca1 = StructureTools.getAtomCAArray(struc_ini);
			Atom[] ca2 = StructureTools.getAtomCAArray(struc_fin);
			
			// get default parameters
			FatCatParameters params = new FatCatParameters();
			
			
			AFPChain afpChain = algorithm.align(ca1,ca2,params);            
			
			//afpChain.setName1(name1);
			//afpChain.setName2(name2);

			currentRMSd = afpChain.getChainRmsd();
	} catch (Exception e) {
		e.printStackTrace();
	}
	}

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
	PDBfromFASTA.writePDB(dir + name + "_" + args[1] + ".pdb", struc_fin);
}

}
