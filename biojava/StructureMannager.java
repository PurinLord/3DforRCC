package rccto3d.optimisation;

import java.util.Vector;
import java.util.Random;
import java.lang.Math;
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
import org.biojava.bio.structure.align.ce.CeParameters;
import org.biojava.bio.structure.align.ce.CeMain;
import org.biojava.bio.structure.align.seq.SmithWaterman3Daligner;

import rccto3d.*;
import rccto3d.Trans;
import rccto3d.PDBfromFASTA;

//import org.SmithWaterman3Daligner;

public class StructureMannager {

	private int initialStructureType;
	// 0 fasta
	// 1 pdb
	// 2 get target
	private int targetStructureType;
	// 0 un target
	// 1 multi targat
	private int energyFunction;
	// 0 SmithWaterman
	// 1 SPalignNS
	// 2 RCC
	private int modifyMethod;
	// 0 min-max change
	// 1 target mold
	// 2 all random
	// 3 substitutor

	//factor de cambio para Phi y Psi
	private double initPhi;
	private double initPsi;
	private double minPhi;
	private double minPsi;
	private double cambioPhi;
	private double cambioPsi;
	private double angleSteps = 0;

	private String fastaID = "a";
	private String fileDir = "";

	private Structure struc_ini = null;
	private Structure[] struc_target;
	private Structure struc_fit = null;
	private Structure struc_model = null;
	private Structure struc_best = null;

	private long initSeed;

	private int multiEnergy = 1;

	private Substitutor sub = null;

	//TODO no se usa para el Calculo de energía
	private int targetRCC[][]; 

	private Structure target = null;
	
	public StructureMannager(int initialStructureType, int targetStructureType, int energyFunction, int modifyMethod){
		this.initialStructureType = initialStructureType;
		this.targetStructureType = targetStructureType;
		this.modifyMethod = modifyMethod;
		this.energyFunction = energyFunction;
	}

	public void setConditions(double anguloInicial, double anguloFinal, String fileDir, long initSeed){

		this.initPhi = this.cambioPhi = this.initPsi = this.cambioPsi = anguloInicial;
		this.minPhi = this.minPsi = anguloFinal;

		this.fileDir = fileDir;
		this.initSeed = initSeed;
	}

	public void setFastaID(String id){
		fastaID = id;
	}

	public void setSubsitutor(Substitutor sub){
		this.sub = sub;
	}

	public String loadStructures(String[] fileName, int verbos){
		// Load Protein Sequence
		PDBfromFASTA pff = null;
		String pdb = null;
		String outString = null;
		
		//load initial
		switch (initialStructureType){
			case 0: //fasta
				try{
					pff = new PDBfromFASTA();
					pdb = pff.pdbFromFile(fileName[0], fastaID);
				
					// Build initial Protein 3D structure
					struc_ini = pff.randomShapeProtein(pdb, initSeed);
					
				}catch(Exception e){
					e.printStackTrace();
				}
				break;

			case 1: //pdb
				struc_ini = PDBfromFASTA.readPDB(fileName[0]);
				break;

			case 2: //no init Copy target
				if(targetStructureType == 0){
					struc_ini = PDBfromFASTA.readPDB(fileName[0]);
					struc_ini = alterConformationAll(struc_ini, 180, 180);
					target = PDBfromFASTA.readPDB(fileName[0]);
				}
				// TODO elegir un objetibo al azar
				break;
		}

		//Load target
		switch (targetStructureType){
			case 0:
				struc_target = new Structure[1];
				if(initialStructureType != 2){
					struc_target[0] = PDBfromFASTA.readPDB(fileName[1]);
				}else{
					struc_target[0] = PDBfromFASTA.readPDB(fileName[0]);
				}
				break;
			case 1:
				int offset = 0;
				if(initialStructureType != 2){
					offset = 1;
				}
				struc_target = new Structure[fileName.length -offset];
				for(int i = 0;i<struc_target.length;i++){
					struc_target[i] = PDBfromFASTA.readPDB(fileName[i+offset]);
				}
				break;
		}

		//substitutor
		if(modifyMethod == 3){
			if(sub == null){
				sub = new Substitutor(struc_target[0]);
				// TODO elegir un objetibo al azar
				sub.randomInitialize();
			}
			if(initialStructureType != 1){
				struc_ini = sub.fakeSubstitute(struc_ini);
			}
		}
		struc_best = (Structure)struc_ini.clone();
		PDBfromFASTA.writePDB(fileDir + "struc_ini.pdb", struc_ini);

		if(targetStructureType == 0){
			PDBfromFASTA.writePDB(fileDir + "struc_target.pdb", struc_target[0]);
		}
		if(targetStructureType == 1){
			int i = 0;
			for(Structure struc : struc_target){
				PDBfromFASTA.writePDB(fileDir + "struc_target_"+i+".pdb", struc);
				i ++;
			}
		}

		if(energyFunction == 2 || multiEnergy > 1){
			targetRCC = new int[struc_target.length][26];
			int i = 0;
			for(Structure struc : struc_target){
				targetRCC[i] = calcRCC(struc);
				i++;
			}
		}

		if (verbos > 0){
			outString = fileName[0] + " " + fileName[1]; //TODO un for para poner todo filename
			outString += "\nangulo Inicial " + initPhi + " angulo final " + minPhi + " initSeed " + initSeed;
			if(modifyMethod == 3){
				Vector<Vector<Integer>> divition = sub.getDivition();
				outString += "\n" + divition;
			}
			if (verbos > 1){
				if(energyFunction == 2 || multiEnergy > 1){
					for(int[] vect : targetRCC){
						for(int i : vect){System.out.print(i + " ");}
					}
				}
				System.out.println(outString);
			}
			return outString;
		}
		return null;
	}

	public static int[] calcRCC(String s1, String tmpdir){
		int rcc[] = new int[26];
		try{
		String s2 = "A";
		Process p = Runtime.getRuntime().exec("python create_26dvRCC.py "+s1+" "+s2 +" "+tmpdir);
		BufferedReader in = new BufferedReader(new InputStreamReader(p.getInputStream()));
		String ret = in.readLine();
		String val = "";
		int j = 0;
		for(int i = 0; i < ret.length(); i++){
			char c = ret.charAt(i);
			if(c==','){
				rcc[j] = Integer.valueOf(val);
				val = "";
				j++;
			}else{
				val+=c;
			}
		}
		rcc[j] = Integer.valueOf(val);
		}catch(Exception e){
  		e.printStackTrace();
		}
		return rcc;
	}

	public int[] calcRCC(String s1){
		return calcRCC(s1, fileDir);
	}

	//TODO
	public int[] calcRCC(Structure struc) {
		int[] result = new int[26];
		int i=0;
		
		for(i=0; i<26; i++) result[i]=0;
		
		return result;
	} // end calcRCC

	public double calcRMSD(Structure struc_target, Structure struc_current){
		double currentRMSd = 0.0;
		try {
			StructureAlignment algorithm  = StructureAlignmentFactory.getAlgorithm(FatCatRigid.algorithmName);

			Atom[] ca1 = StructureTools.getAtomCAArray(struc_target);
			Atom[] ca2 = StructureTools.getAtomCAArray(struc_current);
			
			FatCatParameters params = new FatCatParameters();
			
			AFPChain afpChain = algorithm.align(ca1,ca2,params);            
			
			currentRMSd = afpChain.getChainRmsd();
		} catch (Exception e) {
		e.printStackTrace();
		}
		return currentRMSd;
	}

	public double calcBackboneRMSD(Structure struc_target, Structure struc_current){
		double currentRMSd = 0.0;
		try {
			StructureAlignment algorithm  = StructureAlignmentFactory.getAlgorithm(FatCatRigid.algorithmName);

			Atom[] ca1 = StructureTools.getBackboneAtomArray(struc_target);
			Atom[] ca2 = StructureTools.getBackboneAtomArray(struc_current);
			
			FatCatParameters params = new FatCatParameters();
			
			AFPChain afpChain = algorithm.align(ca1,ca2,params);            
			
			currentRMSd = afpChain.getChainRmsd();
		} catch (Exception e) {
		e.printStackTrace();
		}
		return currentRMSd;
	}

	public double calcRMSD2(Structure struc_target, Structure struc_current){
		double currentRMSd = 0.0;
		try {
 			StructureAlignment algorithm  = StructureAlignmentFactory.getAlgorithm(CeMain.algorithmName);
 			
			Atom[] ca1 = StructureTools.getAtomCAArray(struc_target);
			Atom[] ca2 = StructureTools.getAtomCAArray(struc_current);
 			
 			// get default parameters
 			CeParameters params = new CeParameters();
 			// set the maximum gap size to unlimited 
 			//params.setMaxGapSize(-1);
 			
 			AFPChain afpChain = algorithm.align(ca1,ca2,params);            

			currentRMSd = afpChain.getChainRmsd();
		}catch (Exception e) {
			e.printStackTrace();
		}
		return currentRMSd;
	}

	public double calcBackboneRMSD2(Structure struc_target, Structure struc_current){
		double currentRMSd = 0.0;
		try {
 			StructureAlignment algorithm  = StructureAlignmentFactory.getAlgorithm(CeMain.algorithmName);
 			
			Atom[] ca1 = StructureTools.getBackboneAtomArray(struc_target);
			Atom[] ca2 = StructureTools.getBackboneAtomArray(struc_current);
 			
 			// get default parameters
 			CeParameters params = new CeParameters();
 			// set the maximum gap size to unlimited 
 			params.setMaxGapSize(-1);
 			
 			AFPChain afpChain = algorithm.align(ca1,ca2,params);            

			currentRMSd = afpChain.getChainRmsd();
		}catch (Exception e) {
			e.printStackTrace();
		}
		return currentRMSd;
	}

	public double calcRMSD_SmithWaterman(Structure struc_current, Structure struc_target){
		double currentRMSd = 0.0;
		try {
			Atom[] ca1 = StructureTools.getAtomCAArray(struc_current);
			Atom[] ca2 = StructureTools.getAtomCAArray(struc_target);
			SmithWaterman3Daligner swAligner = new SmithWaterman3Daligner();
 			
 			AFPChain afpChain = swAligner.align(ca1,ca2);            

			currentRMSd = afpChain.getChainRmsd();
		}catch (Exception e) {
			e.printStackTrace();
		}
		return currentRMSd*100;
	}

	public double calcRMSD_SPalignNS(String dir1, String dir2){
		double rmsd=0;
		try{
		String s1 = "-rmsOnly -pair";
		Process p = Runtime.getRuntime().exec("./SPalignNS "+s1+" "+dir1+" "+dir2);
		BufferedReader in = new BufferedReader(new InputStreamReader(p.getInputStream()));
		rmsd = Double.parseDouble(in.readLine());
		}catch(Exception e){
  		e.printStackTrace();
		}
		return rmsd;
	}

	public double calcRMSD_SPalignNS(Structure struc_current, Structure struc_target){
		PDBfromFASTA.writePDB(fileDir + "struc_1.pdb", struc_current);
		PDBfromFASTA.writePDB(fileDir + "struc_2.pdb", struc_target);
		return calcRMSD_SPalignNS(fileDir + "struc_1.pdb", fileDir + "struc_2.pdb");
	}

 	//static {
  	//System.loadLibrary("SPalig");
	//}
	//public native double calcRMSD_SPalignNS_fast(String pdb_struc_current, String pdb_struc_target);
	//public double calcRMSD_SPalignNS_fast(Structure struc_current, Structure struc_target){
	//	//if(target_pdb == null){
	//	//	target_pdb = struc_target.toPDB();
	//	//}
	//	String s1 = struc_target.toPDB();
	//	String s2 = struc_current.toPDB();
	//	return 100*calcRMSD_SPalignNS_fast(s1, s2);
	//}

	public double calcRCCSimilarity(int[] rcc1, int[] rcc2){
		int sum = 0;
		int i=0;
		for(i=0; i<26; i++){
			sum += Math.pow(rcc1[i] - rcc2[i], 2.0);
		}
		return Math.sqrt(sum);
	}

	public double calcEnergy(Structure[] target, Structure model){
		double[] energy = new double[target.length];
		int i = 0;
		for(Structure struc : target){
			switch(energyFunction){
				case 0:
					energy[i] = calcRMSD_SmithWaterman(struc, model);
					break;

				case 1:
					energy[i] = calcRMSD_SPalignNS(struc, model);
					break;

				case 2:
					PDBfromFASTA.writePDB(fileDir + "struc_1.pdb", struc);
					int RCC1[] = calcRCC(fileDir + "struc_1.pdb");
					PDBfromFASTA.writePDB(fileDir + "struc_2.pdb", model);
					int RCC2[] = calcRCC(fileDir + "struc_2.pdb");
					energy[i] = calcRCCSimilarity(RCC1, RCC2);;
					break;
			}
			i++;
		}
		if(targetStructureType != 1){
			return energy[0];
		}else{
			double mean = 0;
			mean = mean/energy.length;
			return mean;
		}
	}

	public double multiCalcEnergy(double[] energy){
		return 0;
	}

	public double calcEnergy(){
		//initialize fit for initial energy
		if(struc_fit == null){
			struc_fit = struc_ini;
			struc_best = (Structure)struc_fit.clone();
			return calcEnergy(struc_target, struc_ini);
		}
		return calcEnergy(struc_target, struc_model);
	}

	public void alterConformation(){
		switch(modifyMethod){
			case 0:
				//alterConformationSome(minCambio,maxCambio, struc_fit);
				break;

			case 1:
				struc_model = alterConformationAll(struc_fit, target);
				break;

			case 2:
				struc_model = alterConformationAll(struc_fit);
				break;

			case 3:
				struc_model = alterConformationParts(struc_fit);
				break;
		}
	}

	public void acceptModel(){
		struc_fit = (Structure)struc_model.clone();
	}

	public void setFitAsBeast(){
		struc_best = (Structure)struc_model.clone();
	}

	public void saveBest(){
		double best = calcEnergy(struc_target, struc_best);
		PDBfromFASTA.writePDB(fileDir + "sol_" + best + ".pdb", struc_best);
	}

	//TODO no se que pedo con esto o.O
	public void printEnergy(){
		//if(energyFunction == 2){
		//	for(int i : modelRCC){System.out.print(i + " ");}
		//	System.out.print("\n");
		//}
		
		System.out.println(calcEnergy(struc_target, struc_model));
		if(multiEnergy > 1){
			double energy = calcEnergy(struc_target, struc_model);
			System.out.println(energy);
		}
	}

	public Structure alterConformationSome(int minCambio, int maxCambio, Structure struc){
		Random rdm =  new Random(System.currentTimeMillis());
		Structure alterStruc = (Structure)struc.clone();
		Chain chain = alterStruc.getChain(0);
		int largo = chain.getAtomLength();
		int to = minCambio + (int)(rdm.nextDouble()*(maxCambio - minCambio));
		for(int i=0; i<to; i++){
			int dAmino = (int) (largo * rdm.nextDouble());
			double dPsi = Math.round(cambioPsi * 2*(rdm.nextDouble() - 0.5));
			double dPhi = Math.round(cambioPhi * 2*(rdm.nextDouble() - 0.5));
			Trans.rotatePsi(chain, dAmino, dPsi);
			Trans.rotatePhi(chain, dAmino, dPhi);
		}
		return alterStruc;
	}

	//TODO ROUND
	public Structure alterConformationAll(Structure struc){
		Random rdm =  new Random(System.currentTimeMillis());
		Structure alterStruc = (Structure)struc.clone();
		Chain chain = alterStruc.getChain(0);
		int largo = chain.getAtomLength();
		for(int i=0; i<largo; i++){
			double dPsi = cambioPsi * 2*(rdm.nextDouble() - 0.5);
			double dPhi = cambioPhi * 2*(rdm.nextDouble() - 0.5);
			//System.out.println(dPhi +" "+ dPsi);
			Trans.rotatePsi(chain, i, dPsi);
			Trans.rotatePhi(chain, i, dPhi);
		}
		return alterStruc;
	}

	public Structure alterConformationParts(Structure struc){
		Random rdm =  new Random(System.currentTimeMillis());
		Structure alterStruc = (Structure)struc.clone();
		Chain chain = alterStruc.getChain(0);
		int largo = chain.getAtomLength();
		double dPsi;
		double dPhi;
		int index = 0;
		Vector<Vector<Integer>> divition = sub.getDivition();
		Vector<Integer> segment = divition.elementAt(index);
		for(int i=0; i<largo; i++){
			if(i == segment.elementAt(0)){
				i += segment.elementAt(1)+1;
				index++;
				if(index < divition.size()){
					segment = divition.elementAt(index);
				}else{
					break;
				}
			}
			dPsi = cambioPsi * 2*(rdm.nextDouble() - 0.5);
			dPhi = cambioPhi * 2*(rdm.nextDouble() - 0.5);
			//System.out.println(dPhi +" "+ dPsi);
			Trans.rotatePsi(chain, i, dPsi);
			Trans.rotatePhi(chain, i, dPhi);
		}
		return alterStruc;
	}

	public Structure alterConformationAll(Structure struc, double deltaPsi, double deltaPhi){
		Random rdm =  new Random(System.currentTimeMillis());
		Structure alterStruc = (Structure)struc.clone();
		Chain chain = alterStruc.getChain(0);
		int largo = chain.getAtomLength();
		for(int i=0; i<largo; i++){
			double dPsi = deltaPsi * 2*(rdm.nextDouble() - 0.5);
			double dPhi = deltaPhi * 2*(rdm.nextDouble() - 0.5);
			//System.out.println(i +" "+ dPhi +" "+ dPsi);
			Trans.rotatePsi(chain, i, dPsi);
			Trans.rotatePhi(chain, i, dPhi);
		}
		return alterStruc;
	}

	public Structure alterConformationAll(Structure struc, Structure target){
		Random rdm =  new Random(System.currentTimeMillis());
		Chain chainT = target.getChain(0);
		Structure alterStruc = (Structure)struc.clone();
		Chain chain = alterStruc.getChain(0);
		int largo = chain.getAtomLength();
		AminoAcid a1;
		AminoAcid a2;
		for(int i=0; i <largo; i++){
			int dAmino = (int) ((largo -1) * rdm.nextDouble());
			try{
			a1 = (AminoAcid)chainT.getAtomGroup(dAmino);
			a2 = (AminoAcid)chainT.getAtomGroup(dAmino + 1);
			double dPsi = Trans.getPhi(a1, a2);
			double dPhi = Trans.getPsi(a1, a2);
			Trans.rotatePsi(chain, i, dPsi);
			Trans.rotatePhi(chain, i, dPhi);
    	}catch (Exception e) {
				//System.out.println(">" + largo +" "+ i);
    		e.printStackTrace();
    	}
		}
		return alterStruc;
	}

	//
	//Random de 0 error
	//al dividir salen 0 donde no deberían
	//todo mayor que 0!!!!!

	public void angleCooling(){
		cambioPhi = initPhi*Math.exp(-angleSteps/(-1.0/Math.log(minPhi/initPhi)));
		cambioPsi = cambioPhi;
	}
	public void addAngleStep(int searchStepsCicle, int searchStepsTotal){
		angleSteps += (double)searchStepsCicle/searchStepsTotal;
	}

	public double tempCalc(){
		//Structure alterStruc = sm.alterConformationAll(alterStruc, 180, 180);
		Structure alterStruc = alterConformationAll(struc_ini, 180, 180);
		return calcEnergy(struc_target,alterStruc);
	}

	public static void main(String[] args){
	}
}
