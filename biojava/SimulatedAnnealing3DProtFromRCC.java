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
	double temp;
	int searchSteps;

  // Cooling rate
  double coolingRate;

	//factor de cambio para Phi y Psi
	double cambioPhi;
	double cambioPsi;

	private String dirStartStruc;
	private String fastaID = "a";
	private String dirTargetStruc;
	private String fileDir = "";

	private long initSeed;

	private String startPdb;
	private int targetRCC[]; 

	public SimulatedAnnealing3DProtFromRCC(String dirStartStruc, String dirTargetStruc){
		temp = 10000;
		searchSteps = 3;

  	coolingRate = 0.003;

		cambioPhi = 36.0;
		cambioPsi = 36.0;
		this.dirStartStruc = dirStartStruc;
		this.dirTargetStruc = dirTargetStruc;
		initSeed = (long)Math.random();
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

	public SimulatedAnnealing3DProtFromRCC(double temp, int searchSteps, double coolingRate, double cambioPhi, double cambioPsi,
																					String startPdb, int targetRCC[], String fileDir, long initSeed){
		this.temp = temp;
		this.searchSteps = searchSteps;

  	this.coolingRate = coolingRate;

		this.cambioPhi = cambioPhi;
		this.cambioPsi = cambioPsi;

		this.startPdb = startPdb;
		this.targetRCC = targetRCC;
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

	public double calcSimilarity(int[] rcc1, int[] rcc2){
		int sum = 0;
    int i=0;
		for(i=0; i<26; i++){
			sum += Math.pow(rcc1[i] - rcc2[i], 2.0);
		}
		return Math.sqrt(sum);
  }// end calcSimilarity

	public Structure alterConformation(int minCambio, int maxCambio, Structure struc){
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

	public Structure alterConformationAll(Structure struc){
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

	static public Structure alterConformationAll(Structure struc, double deltaPsi, double deltaPhi){
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

	public Structure alterConformationAll(Structure struc, Structure target){
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
	public float run(int verbos){
		//Mesure time
		long elapsedTime = 0;
		if (verbos > 0){
			elapsedTime = System.nanoTime();
			System.out.println(dirStartStruc + " " + dirTargetStruc + 
				"\ntemp " + temp + " coolRate " + coolingRate + " searchSteps " + searchSteps +
				"\ncamio Phi " + cambioPhi + " cambio Psi " + cambioPsi);// + " min cambio " + minCambio + " max cambio " +  maxCambio);
		}
		// Load Protein Sequence
		PDBfromFASTA pff = null;
		String pdb = null;
		Structure struc_ini = null;
		
		Structure struc_model = null;
		Structure target = null;
		double currentEnergy = 0, neighbourEnergy = 0;
		double distance_ini = 0;
		
		int currentRCC[]; 
		
		if(dirStartStruc != null){
			if(dirStartStruc.length() > 3){
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
		}else{
			try{
				pff = new PDBfromFASTA();
				struc_ini = pff.randomShapeProtein(startPdb, initSeed);
			}catch(Exception e){
				e.printStackTrace();
			}
		}

		
		// Declaration of variables to store energy values
		if(targetRCC != null){
			targetRCC = calcRCC(dirTargetStruc, fileDir);
		}
		if (verbos > 0){
			for(int i : targetRCC){System.out.print(i + " ");}
			System.out.print("\n");
		}
		
		// Initialize intial solution
		distance_ini = calcSimilarity(targetRCC, calcRCC(fileDir + "struc_ini.pdb", fileDir));
		
		if (verbos > 0){
			System.out.println("Initial solution distance: " + distance_ini);
		}
		double distance_neighbor=0.0;
		
		// Set as current best
		double best = distance_ini;
		    
		// Create new neighbour 3d model
		struc_ini = PDBfromFASTA.readPDB(fileDir + "struc_ini.pdb");
		 
		// Loop until system has cooled
		for (;temp > 1;temp = stdCooling(temp)) 
		{
			for(int step = 0; step < searchSteps; step++){
				currentRCC = calcRCC(fileDir + "struc_ini.pdb", fileDir);
				currentEnergy = calcSimilarity(targetRCC, currentRCC);;
				// Get a random conformation for this new neighbor
				if(target != null){
					struc_model = alterConformationAll(struc_ini, target);
				}else{
					struc_model = alterConformationAll(struc_ini);
				}
				PDBfromFASTA.writePDB(fileDir + "struc_model.pdb", struc_model);
							
				currentRCC = calcRCC(fileDir + "struc_model.pdb", fileDir);
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
					PDBfromFASTA.writePDB(fileDir + "sol_" + best + ".pdb", struc_ini);
					if (verbos > 1){
						if (verbos > 0){
							System.out.println(temp + "> " + best);
						}
							System.out.println(temp + "> " + best + " time: " + elapsedTime/3600000000000.0);
					}
					if(best == 0){
						temp = 0;
					}
				}
			
				if (verbos > 1){
					for(int i : currentRCC){System.out.print(i + " ");}
					System.out.print("\n");
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
		simA.run(2);
	}
}

