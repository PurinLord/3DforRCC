package rccto3d.optimisation;

import java.lang.Math;
import java.io.*;

import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.AminoAcid;

import rccto3d.*;
import rccto3d.Trans;
import rccto3d.PDBfromFASTA;

public class SimulatedAnnealing3DProtFromRCC 
{

	// Set initial temp
	static double temp = 10000;
	static int searchSteps = 300;

  // Cooling rate
  static double coolingRate = 0.003;

	//factor de cambio para Phi y Psi
	static double cambioPhi = 36.0;
	static double cambioPsi = 36.0;

	//Candidad de ángulos phi psi que cambian
	static int minCambio = 1;
	static int maxCambio = 70;

	//Calculate the acceptance probability
	public static double acceptanceProbability(double energy, double newEnergy, double temperature) 
	{
		// If the new solution is better, accept it
		if (newEnergy < energy) return 1.0;
		// If the new solution is worse, calculate an acceptance probability
		return Math.exp((energy - newEnergy) / temperature);
	}

	public static int[] calcRCC(String s1){
		int rcc[] = new int[26];
		try{
		String s2 = "A";
		Process p = Runtime.getRuntime().exec("python create_26dvRCC.py "+s1+" "+s2);
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

	public static int[] calcRCC(Structure struc)
	{
		int[] result = new int[26];
		int i=0;
		
		for(i=0; i<26; i++) result[i]=0;
		
		return result;
	} // end calcRCC

	public static double calcSimilarity(int[] rcc1, int[] rcc2)
	{
		int sum = 0;
    int i=0;
		for(i=0; i<26; i++){
			sum += Math.pow(rcc1[i] - rcc2[i], 2.0);
		}
		return Math.sqrt(sum);
  }// end calcSimilarity

	public static Structure alterConformation(Structure struc){
		Chain chain = struc.getChain(0);
		int largo = chain.getAtomLength();
		int to = minCambio + (int)(Math.random()*(maxCambio - minCambio));
		for(int i=0; i<to; i++){
			int dAmino = (int) (largo * Math.random());
			double dPsi = Math.round(cambioPsi * 2*(Math.random() - 0.5));
			double dPhi = Math.round(cambioPhi * 2*(Math.random() - 0.5));
			Trans.rotatePsi(chain, dAmino, dPsi);
			Trans.rotatePhi(chain, dAmino, dPhi);
		}
		return struc;
	}

	public static Structure alterConformationAll(Structure struc){
		Chain chain = struc.getChain(0);
		int largo = chain.getAtomLength();
		for(int i=0; i<largo; i++){
			double dPsi = cambioPsi * 2*(Math.random() - 0.5);
			double dPhi = cambioPhi * 2*(Math.random() - 0.5);
			Trans.rotatePsi(chain, i, dPsi);
			Trans.rotatePhi(chain, i, dPhi);
		}
		return struc;
	}

	public static Structure alterConformationAll(Structure struc, double deltaPsi, double deltaPhi){
		Chain chain = struc.getChain(0);
		int largo = chain.getAtomLength();
		for(int i=0; i<largo; i++){
			double dPsi = deltaPsi * 2*(Math.random() - 0.5);
			double dPhi = deltaPhi * 2*(Math.random() - 0.5);
			Trans.rotatePsi(chain, i, dPsi);
			Trans.rotatePhi(chain, i, dPhi);
		}
		return struc;
	}

	public static Structure alterConformationAll(Structure struc, Structure target){
		Chain chainT = target.getChain(0);
		Chain chain = struc.getChain(0);
		int largo = chain.getAtomLength();
		AminoAcid a1;
		AminoAcid a2;
		for(int i=0; i <largo; i++){
			int dAmino = (int) ((largo -1) * Math.random());
			try{
			a1 = (AminoAcid)chain.getAtomGroup(dAmino);
			a2 = (AminoAcid)chain.getAtomGroup(dAmino + 1);
			double dPsi = Trans.getPhi(a1, a2);
			double dPhi = Trans.getPsi(a1, a2);
			Trans.rotatePsi(chain, i, dPsi);
			Trans.rotatePhi(chain, i, dPhi);
    	}catch (Exception e) {
				System.out.println(">" + largo +" "+ i);
    		e.printStackTrace();
    	}
		}
		return struc;
	}

	public static double stdCooling(double temp){
		return temp * 1-coolingRate;
	}

	public static double linearCooling(double temp){
		return temp - coolingRate;
	}

	public static double slowCooling(double temp){
		return temp/Math.log(temp);
	}

	public static void main(String[] args)
	{
		//Mesure time
		long elapsedTime = System.nanoTime();
		System.out.println(args[0] + " " + args[1] + " " + args[2] + 
				"\ntemp " + temp + " coolRate " + coolingRate + " searchSteps " + searchSteps +
				"\ncamio Phi " + cambioPhi + " cambio Psi " + cambioPsi + " min cambio " + minCambio + " max cambio " +  maxCambio);
		// Load Protein Sequence
		PDBfromFASTA pff = null;
		String fileDir = "";
		if(args.length > 3){
			pff = new PDBfromFASTA(args[3]);
			fileDir = args[3];
		}else{
			pff = new PDBfromFASTA();
			fileDir = "out/";
		}
		String pdb = null;
		Structure struc_ini = null;
		
		Structure struc_model = null;
		Structure target = null;
		double currentEnergy = 0, neighbourEnergy = 0;
		double distance_ini = 0;
		
		int targetRCC[]; 
		int currentRCC[]; 
		
		if(args[0].length() > 3){
			if(args[0].substring(args[0].length() - 3).equals(".fa")){
				try{
					pdb = pff.pdbFromFile(args[0], args[1]);
					//ProteinSequence seq = new ProteinsSequence(args[0]);
				
					// Build initial Protein 3D structure
					struc_ini = pff.makeProtein(pdb);
				}catch(Exception e){
					e.printStackTrace();
				}
				PDBfromFASTA.writePDB(fileDir + "struc_ini.pdb", struc_ini);
			}
			if(args[0].substring(args[0].length() - 4).equals(".pdb")){
				struc_ini = PDBfromFASTA.readPDB(args[0]);
				PDBfromFASTA.writePDB(fileDir + "struc_ini.pdb", struc_ini);
			}
		}else{
			struc_ini = PDBfromFASTA.readPDB(args[2]);
			target = PDBfromFASTA.readPDB(args[2]);
			PDBfromFASTA.writePDB(fileDir + "struc_ini.pdb", struc_ini);
		}
		
		
		// Declaration of variables to store energy values
		targetRCC = calcRCC(args[2]);
		for(int i : targetRCC){System.out.print(i + " ");}
		System.out.print("\n");
		
		// Initialize intial solution
		distance_ini = calcSimilarity(targetRCC, calcRCC(fileDir + "struc_ini.pdb"));
		
		System.out.println("Initial solution distance: " + distance_ini);
		double distance_neighbor=0.0;
		
		// Set as current best
		double best = distance_ini;
		    
		// Create new neighbour 3d model
		struc_ini = PDBfromFASTA.readPDB(fileDir + "struc_ini.pdb");
		 
		// Loop until system has cooled
		for (;temp > 1;temp = slowCooling(temp)) 
		{
		 
			for(int step = 0; step < searchSteps; step++){
				currentRCC = calcRCC(fileDir + "struc_ini.pdb");
				currentEnergy = calcSimilarity(targetRCC, currentRCC);;
				// Get a random conformation for this new neighbor
				if(target != null){
					struc_model = alterConformationAll(struc_ini, target);
				}else{
					struc_model = alterConformationAll(struc_ini);
				}
				PDBfromFASTA.writePDB(fileDir + "struc_model.pdb", struc_model);
							
				currentRCC = calcRCC(fileDir + "struc_model.pdb");
				// Get energy of solution
				neighbourEnergy = calcSimilarity(targetRCC, currentRCC);;
			
				// Decide if we should accept the neighbour
				if (acceptanceProbability(currentEnergy, neighbourEnergy, temp) > Math.random()) {
					struc_ini = struc_model;
					PDBfromFASTA.writePDB(fileDir + "struc_ini.pdb", struc_ini);
				}
			
				// Keep track of the best solution found
				if (neighbourEnergy < best) 
				{
					best = neighbourEnergy;
					PDBfromFASTA.writePDB(fileDir + "best_" + best + ".pdb", struc_ini);
					System.out.println(temp + ": " + best + " time: " + elapsedTime/3600000000000.0);
					if(best == 0){
						temp = 0;
					}
				}
			
				for(int i : currentRCC){System.out.print(i + " ");}
				System.out.print("\n");
			}
		}
		
		System.out.println("Final solution distance: " + best);//.getDistance());
		elapsedTime = System.nanoTime() - elapsedTime;
		System.out.println("Total execution time: " + elapsedTime/3600000000000.0);
  }
}

