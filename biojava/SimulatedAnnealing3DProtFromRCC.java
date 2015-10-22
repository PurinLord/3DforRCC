package rccto3d.optimisation;

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

import rccto3d.*;
import rccto3d.Trans;
import rccto3d.PDBfromFASTA;

public class SimulatedAnnealing3DProtFromRCC 
{

	// Set initial temp
	private double temp;
	private int searchSteps;

  // Cooling rate
  private double coolingRate;

	//factor de cambio para Phi y Psi
	private double cambioPhi;
	private double cambioPsi;

	private String dirStartStruc;
	private String fastaID = "a";
	private String dirTargetStruc;
	private String fileDir = "";

	private long initSeed;

	private int targetRCC[]; 

	Structure target = null;
	//steps = -1/log((1/temp), 1+coolingRate)

	public SimulatedAnnealing3DProtFromRCC(String dirStartStruc, String dirTargetStruc){
		temp = 10000;
		searchSteps = 1;
		//searchSteps = 300;

  	coolingRate = 0.003;
  	//coolingRate = 0.60189;

		cambioPhi = 0.04;
		cambioPsi = 0.04;
		this.dirStartStruc = dirStartStruc;
		this.dirTargetStruc = dirTargetStruc;
		fileDir = "out/";
		initSeed = (long)(1000*Math.random());
	}

	public SimulatedAnnealing3DProtFromRCC(double temp, int searchSteps, double coolingRate, double cambioPhi, double cambioPsi,
																					String dirStartStruc, String dirTargetStruc, String fileDir, long initSeed){
		this.temp = temp;
		this.searchSteps = searchSteps;

  	this.coolingRate = coolingRate;

		this.cambioPhi = cambioPhi;
		this.cambioPsi = cambioPsi;

		this.dirStartStruc = dirStartStruc;
		this.dirTargetStruc = dirTargetStruc;
		this.fileDir = fileDir;
		this.initSeed = initSeed;
	}

	public void setFastaID(String id){
		fastaID = id;
	}

	//Calculate the acceptance probability
	public double acceptanceProbability(double energy, double newEnergy, double temperature) {
		// If the new solution is better, accept it
		if (newEnergy < energy) return 1.0;
		// If the new solution is worse, calculate an acceptance probability
		return Math.exp((energy - newEnergy) / temperature);
	}

	public double rawAcceptanceProbability(double energy, double newEnergy, double blhe) {
		if (newEnergy < energy) return 1.0;
		return 0.0;
	}

	public static int[] calcRCC(String s1){
		return calcRCC(s1, "");
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

	//TODO
	public int[] calcRCC(Structure struc)
	{
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

	public double calcRMSD2(Structure struc_target, Structure struc_current){
		double currentRMSd = 0.0;
		try {
 			StructureAlignment algorithm  = StructureAlignmentFactory.getAlgorithm(CeMain.algorithmName);
 			
			Atom[] ca1 = StructureTools.getAtomCAArray(struc_target);
			Atom[] ca2 = StructureTools.getAtomCAArray(struc_current);
 			
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

	public double calcSimilarity(int[] rcc1, int[] rcc2){
		int sum = 0;
    int i=0;
		for(i=0; i<26; i++){
			sum += Math.pow(rcc1[i] - rcc2[i], 2.0);
		}
		return Math.sqrt(sum);
  }// end calcSimilarity

	public Structure alterConformation(int minCambio, int maxCambio, Structure struc){
		Structure alterStruc = (Structure)struc.clone();
		Chain chain = alterStruc.getChain(0);
		int largo = chain.getAtomLength();
		int to = minCambio + (int)(Math.random()*(maxCambio - minCambio));
		for(int i=0; i<to; i++){
			int dAmino = (int) (largo * Math.random());
			double dPsi = Math.round(cambioPsi * 2*(Math.random() - 0.5));
			double dPhi = Math.round(cambioPhi * 2*(Math.random() - 0.5));
			Trans.rotatePsi(chain, dAmino, dPsi);
			Trans.rotatePhi(chain, dAmino, dPhi);
		}
		return alterStruc;
	}

	//TODO ROUND
	public Structure alterConformationAll(Structure struc){
		Structure alterStruc = (Structure)struc.clone();
		Chain chain = alterStruc.getChain(0);
		int largo = chain.getAtomLength();
		for(int i=0; i<largo; i++){
			double dPsi = cambioPsi * 2*(Math.random() - 0.5);
			double dPhi = cambioPhi * 2*(Math.random() - 0.5);
			//System.out.println(dPhi +" "+ dPsi);
			Trans.rotatePsi(chain, i, dPsi);
			Trans.rotatePhi(chain, i, dPhi);
		}
		return alterStruc;
	}

	public Structure alterConformationAll(Structure struc, double deltaPsi, double deltaPhi){
		Structure alterStruc = (Structure)struc.clone();
		Chain chain = alterStruc.getChain(0);
		int largo = chain.getAtomLength();
		for(int i=0; i<largo; i++){
			double dPsi = deltaPsi * 2*(Math.random() - 0.5);
			double dPhi = deltaPhi * 2*(Math.random() - 0.5);
			//System.out.println(i +" "+ dPhi +" "+ dPsi);
			Trans.rotatePsi(chain, i, dPsi);
			Trans.rotatePhi(chain, i, dPhi);
		}
		return alterStruc;
	}

	public Structure alterConformationAll(Structure struc, Structure target){
		Chain chainT = target.getChain(0);
		Structure alterStruc = (Structure)struc.clone();
		Chain chain = alterStruc.getChain(0);
		int largo = chain.getAtomLength();
		AminoAcid a1;
		AminoAcid a2;
		for(int i=0; i <largo; i++){
			int dAmino = (int) ((largo -1) * Math.random());
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

	public double stdCooling(double temperature){
		return temperature * (1-coolingRate);
	}

	public double linearCooling(double temp){
		return temp - coolingRate;
	}

	public double slowCooling(double temp){
		return temp/Math.log(temp);
	}

	//verbos 0 - none
	//			 1 - just solutions
	//			 2 - all
	
	public void initialize(int verbos){
		if (verbos > 0){
			System.out.println(dirStartStruc + " " + dirTargetStruc + 
				"\ntemp " + temp + " coolRate " + coolingRate + " searchSteps " + searchSteps +
				"\ncamio Phi " + cambioPhi + " cambio Psi " + cambioPsi + " initSeed " + initSeed);
		}
		// Load Protein Sequence
		PDBfromFASTA pff = null;
		String pdb = null;
		Structure struc_ini = null;

		if(dirStartStruc != null){
			if(dirStartStruc.substring(dirStartStruc.length() - 3).equals(".fa")){
				try{
					pff = new PDBfromFASTA();
					pdb = pff.pdbFromFile(dirStartStruc, fastaID);
					//ProteinSequence seq = new ProteinsSequence(dirStartStruc);
				
					// Build initial Protein 3D structure
					struc_ini = pff.randomShapeProtein(pdb, initSeed);
				}catch(Exception e){
					e.printStackTrace();
				}
				PDBfromFASTA.writePDB(fileDir + "struc_ini.pdb", struc_ini);
			}
			if(dirStartStruc.substring(dirStartStruc.length() - 4).equals(".pdb")){
				struc_ini = PDBfromFASTA.readPDB(dirStartStruc);
				PDBfromFASTA.writePDB(fileDir + "struc_ini.pdb", struc_ini);
			}
		}else{
			struc_ini = PDBfromFASTA.readPDB(dirTargetStruc);
			target = PDBfromFASTA.readPDB(dirTargetStruc);
			PDBfromFASTA.writePDB(fileDir + "struc_ini.pdb", struc_ini);
		}
		
		// Declaration of variables to store energy values
		targetRCC = calcRCC(dirTargetStruc, fileDir);
		if (verbos > 0){
			for(int i : targetRCC){System.out.print(i + " ");}
			System.out.print("\n");
		}
	}

	public float run(int verbos){
		//Mesure time
		long elapsedTime = 0;
		if (verbos > 0){
			elapsedTime = System.nanoTime();
		}
		
		Structure struc_fit = null;
		Structure struc_model = null;
		double currentEnergy[] = {0.0,0.0};
		double neighbourEnergy[] = {0.0,0.0};
		double distance_ini[] = {0.0,0.0};
		double best;

		Structure struc_ini = PDBfromFASTA.readPDB(fileDir + "struc_ini.pdb");
		Structure struc_target = PDBfromFASTA.readPDB(dirTargetStruc);
		
		int currentRCC[]; 
		int modelRCC[]; 
		
		// Initialize intial solution
		currentRCC = calcRCC(fileDir + "struc_ini.pdb", fileDir);
		distance_ini[0] = calcSimilarity(targetRCC, currentRCC);
		distance_ini[1] = calcRMSD(struc_ini,struc_target);
		
		if (verbos > 0){
			System.out.println("Initial solution distance: " + distance_ini[0] + " " + distance_ini[1]);
		}
		
		// Set as current best
		best = distance_ini[1];
		    
		// Create new neighbour 3d model
		struc_fit = PDBfromFASTA.readPDB(fileDir + "struc_ini.pdb");
		PDBfromFASTA.writePDB(fileDir + "struc_fit.pdb", struc_fit);
		 
		currentEnergy[0] = distance_ini[0];
		currentEnergy[1] = distance_ini[1];

		// Loop until system has cooled
		for (;temp > 1;temp = stdCooling(temp)){
			// Search steps
			for(int step = 0; step < searchSteps; step++){
				// Get a random conformation for this new neighbor
				if(target != null){
					struc_model = alterConformationAll(struc_fit, target);
				}else{
					struc_model = alterConformationAll(struc_fit);
				}
				//PDBfromFASTA.writePDB(fileDir + "struc_model.pdb", struc_model);
				//modelRCC = calcRCC(fileDir + "struc_model.pdb", fileDir);

				// Get energy of solution
				//neighbourEnergy[0] = calcSimilarity(targetRCC, modelRCC);;
				neighbourEnergy[1] = calcRMSD(struc_target, struc_model);
			
				// Decide if we should accept the neighbour
				if (acceptanceProbability(currentEnergy[1], neighbourEnergy[1], temp) > Math.random()) {
					struc_fit = struc_model;
					//PDBfromFASTA.writePDB(fileDir + "struc_fit.pdb", struc_fit);
					//currentEnergy[0] = neighbourEnergy[0];
					currentEnergy[1] = neighbourEnergy[1];
					System.out.println("^ " + acceptanceProbability(currentEnergy[1], neighbourEnergy[1], temp));
				}
			
				// Keep track of the best solution found
				if (neighbourEnergy[1] < best){
					best = neighbourEnergy[1];
					PDBfromFASTA.writePDB(fileDir + "sol_" + best + ".pdb", struc_fit);
					if (verbos > 0){
						if (verbos > 1){
							System.out.println(temp + "> " + best + " time: " + elapsedTime/3600000000000.0);
						}else{
							System.out.println(temp + "> " + best);
						}
					}
					if(best == 0){
						temp = 0;
					}
				}
				if (verbos > 1){
					//for(int i : modelRCC){System.out.print(i + " ");}
					//System.out.print("\n");
					System.out.println(calcRMSD(struc_target, struc_fit) +"\t\t"+ neighbourEnergy[1]);
				}
			}
		}
		
		if (verbos > 0){
			System.out.println("Final solution distance: " + best);//.getDistance());
			elapsedTime = System.nanoTime() - elapsedTime;
			System.out.println("Total execution time: " + elapsedTime/3600000000000.0);
		}
		return (float)best;
  }

	public static void main(String[] args){
		SimulatedAnnealing3DProtFromRCC simA = new SimulatedAnnealing3DProtFromRCC(args[0], args[1]);
		simA.initialize(2);
		simA.run(2);
	}
}
