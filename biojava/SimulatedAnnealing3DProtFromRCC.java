package rccto3d.optimisation;

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

import rccto3d.*;
import rccto3d.Trans;
import rccto3d.PDBfromFASTA;

public class SimulatedAnnealing3DProtFromRCC 
{

	// Set initial temp
	private double temp;
	private double tempPercent = 0.15;
	private int searchStepsTotal;
	private int searchStepsCicle;

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
	// 0 - RCC
	// 1 - RMSD
	private int energyType = 1;
	private boolean dualEnergy = false;

	private int targetRCC[]; 

	Structure target = null;
	//steps = -1/log((1/temp), 1+coolingRate)

	public SimulatedAnnealing3DProtFromRCC(String dirStartStruc, String dirTargetStruc){
		searchStepsCicle = 1;
		this.searchStepsTotal = 3000;
		//searchStepsCicle = 300;
  	//coolingRate = 0.60189;
		this.temp = 0;
		this.coolingRate = 0;

		cambioPhi = 9.0;
		cambioPsi = 9.0;
		this.dirStartStruc = dirStartStruc;
		this.dirTargetStruc = dirTargetStruc;
		fileDir = "out/";
		initSeed = (long)(1000*Math.random());
	}

	public SimulatedAnnealing3DProtFromRCC(int searchStepsTotal, int searchStepsCicle, double cambioPhi, double cambioPsi,
																					String dirStartStruc, String dirTargetStruc, String fileDir, long initSeed){
		this.searchStepsTotal = searchStepsTotal;
		this.searchStepsCicle = searchStepsCicle;
		this.temp = 0;
		this.coolingRate = 0;

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

	public double getTemp(){
		return this.temp;
	}

	public void setTemp(double temp){
		this.temp = temp;
	}

	public void setTempPercent(double tempPercent){
		this.tempPercent = tempPercent;
	}

	public void setCoolRate(double coolingRate){
		this.coolingRate = coolingRate;
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

	public int[] calcRCC(String s1){
		return calcRCC(s1, fileDir);
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

	public double[] calcEnergy(Structure s1, Structure s2){
		double energy[] = {0.0,0.0};
		if(energyType == 0 || dualEnergy){
			PDBfromFASTA.writePDB(fileDir + "struc_1.pdb", s1);
			int RCC1[] = calcRCC(fileDir + "struc_1.pdb");
			PDBfromFASTA.writePDB(fileDir + "struc_2.pdb", s2);
			int RCC2[] = calcRCC(fileDir + "struc_2.pdb");
			energy[0] = calcSimilarity(RCC1, RCC2);;
		}
		if(energyType == 1 || dualEnergy){
			energy[1] = calcBackboneRMSD(s1, s2);
		}
		return energy;
	}

	public double[] calcEnergy(Structure s1, Structure s2, int RCC1[]){
		double energy[] = {0.0,0.0};
		if(energyType == 0 || dualEnergy){
			PDBfromFASTA.writePDB(fileDir + "struc_2.pdb", s2);
			int RCC2[] = calcRCC(fileDir + "struc_2.pdb");
			energy[0] = calcSimilarity(RCC1, RCC2);;
		}
		if(energyType == 1 || dualEnergy){
			energy[1] = calcBackboneRMSD(s1, s2);
		}
		return energy;
	}

	public Structure alterConformation(int minCambio, int maxCambio, Structure struc){
		Structure alterStruc = (Structure)struc.clone();
		Chain chain = alterStruc.getChain(0);
		int largo = chain.getAtomLength();
		Random rdm =  new Random(System.currentTimeMillis());
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
		Structure alterStruc = (Structure)struc.clone();
		Chain chain = alterStruc.getChain(0);
		Random rdm =  new Random(System.currentTimeMillis());
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

	public Structure alterConformationAll(Structure struc, double deltaPsi, double deltaPhi){
		Structure alterStruc = (Structure)struc.clone();
		Chain chain = alterStruc.getChain(0);
		Random rdm =  new Random(System.currentTimeMillis());
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
		Chain chainT = target.getChain(0);
		Structure alterStruc = (Structure)struc.clone();
		Chain chain = alterStruc.getChain(0);
		Random rdm =  new Random(System.currentTimeMillis());
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
	
	public double calcInitialTemp(int verbos, Structure struc_seed){
		PDBfromFASTA.writePDB(fileDir + "calcTmp0.pdb", struc_seed);
		int seedRCC[] = calcRCC(fileDir + "calcTmp0.pdb");
		//int alterRCC[];
		Structure alterStruc = (Structure)struc_seed.clone();
		int count = 0;
		int searchStepsCicle = 50;
		//double max = 0;
		double sum = 0;
		double currentEnergy[]={0.0,0.0};
		while(count < searchStepsCicle){
			alterStruc = alterConformationAll(alterStruc, 180, 180);

			currentEnergy = calcEnergy(struc_seed,alterStruc, seedRCC);
			//PDBfromFASTA.writePDB(fileDir + "calcTmp1.pdb", alterStruc);
			//alterRCC = calcRCC(fileDir + "calcTmp1.pdb");
			//currentEnergy = calcSimilarity(seedRCC, alterRCC);;
			//currentEnergy[1] = calcRMSD(struc_seed,alterStruc);

			//max = Math.max(max, currentEnergy[energyType]);
			sum += currentEnergy[energyType];
			if(verbos > 2){
				System.out.println("energy " + currentEnergy[energyType] + " sum " + sum);
			}
			count++;
		}
		double prom = sum/count;
		if(verbos > 1){
			System.out.println(" promE " + prom);
		}
		//double temp = ((3/2)*prom)+(max/2);
		double temp = prom;//2*max - prom;
		this.temp = temp * tempPercent;
		return temp;
	}

	//TODO
	public double calcInitialTemp(int verbos){
		return calcInitialTemp(verbos, null);
	}

	public void calcCoolRateStd(int verbos, int searchStepsTotal, int searchStepsCicle){
		double steps = searchStepsTotal/(double)searchStepsCicle;
		this.coolingRate = 1-(Math.pow((1/this.temp),(1.0/steps)));
	}
	
	public String initialize(int verbos){
		// Load Protein Sequence
		PDBfromFASTA pff = null;
		String pdb = null;
		Structure struc_ini = null;
		Structure struc_target;
		String outString = null;

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
		struc_target = PDBfromFASTA.readPDB(dirTargetStruc);
		PDBfromFASTA.writePDB(fileDir + "struc_target.pdb", struc_target);
		targetRCC = calcRCC(dirTargetStruc);

		if(temp == 0){
			calcInitialTemp(verbos, struc_target);
		}
		if(coolingRate == 0){
			calcCoolRateStd(verbos, searchStepsTotal, searchStepsCicle);
		}

		if (verbos > 0){
			outString = dirStartStruc + " " + dirTargetStruc + 
				"\ntemp "+temp +" coolRate "+coolingRate+" totalSearchSteps "+searchStepsTotal+" searchStepsCicle "+searchStepsCicle +
				"\ncamio Phi " + cambioPhi + " cambio Psi " + cambioPsi + " initSeed " + initSeed;
			if (verbos > 1){
				for(int i : targetRCC){System.out.print(i + " ");}
				System.out.println(outString);
				System.out.print("\n");
			}
			return outString;
		}
		return null;
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
		Structure struc_target = PDBfromFASTA.readPDB(fileDir + "struc_target.pdb");
		
		int currentRCC[]; 
		int modelRCC[]; 

		Random rdm =  new Random(System.currentTimeMillis());
		
		// Initialize intial solution
		currentRCC = calcRCC(fileDir + "struc_ini.pdb");
		distance_ini[0] = calcSimilarity(targetRCC, currentRCC);
		distance_ini[1] = calcBackboneRMSD(struc_ini,struc_target);
		
		if (verbos > 0){
			System.out.println("Initial solution distance: " + distance_ini[0] + " " + distance_ini[1]);
		}
		
		// Set as current best
		best = distance_ini[energyType];
		    
		// Create new neighbour 3d model
		struc_fit = PDBfromFASTA.readPDB(fileDir + "struc_ini.pdb");
		PDBfromFASTA.writePDB(fileDir + "struc_fit.pdb", struc_fit);
		 
		currentEnergy = distance_ini;

		// Loop until system has cooled
		for (;temp > 1;temp = stdCooling(temp)){
			// Search steps
			for(int step = 0; step < searchStepsCicle; step++){
				// Get a random conformation for this new neighbor
				if(target != null){
					struc_model = alterConformationAll(struc_fit, target);
				}else{
					struc_model = alterConformationAll(struc_fit);
				}
				// Get energy of solution
				neighbourEnergy = calcEnergy(struc_target, struc_model, targetRCC);
				//PDBfromFASTA.writePDB(fileDir + "struc_model.pdb", struc_model);
				//modelRCC = calcRCC(fileDir + "struc_2.pdb");
				//neighbourEnergy[0] = calcSimilarity(targetRCC, modelRCC);
				//neighbourEnergy[1] = calcRMSD(struc_target, struc_model);
			
				// Decide if we should accept the neighbour
				if (acceptanceProbability(currentEnergy[energyType], neighbourEnergy[energyType], temp) > rdm.nextDouble()) {
					struc_fit = struc_model;
					PDBfromFASTA.writePDB(fileDir + "struc_fit.pdb", struc_fit);
					currentEnergy = neighbourEnergy;
					if (verbos > 0){
						if (verbos > 1){
							long time = System.nanoTime() - elapsedTime;
							System.out.println(temp + "> " + currentEnergy[energyType] + " time: " + time/3600000000000.0);
						}else{
							System.out.println(temp + "> " + currentEnergy[energyType]);
						}
					}
				}
			
				// Keep track of the best solution found
				if (neighbourEnergy[energyType] < best){
					best = neighbourEnergy[energyType];
					PDBfromFASTA.writePDB(fileDir + "sol_" + best + ".pdb", struc_fit);
					if(best == 0){
						temp = 0;
					}
				}
				if (verbos > 1){
					if(energyType == 0){
						//for(int i : modelRCC){System.out.print(i + " ");}
						//System.out.print("\n");
					}
					System.out.println(neighbourEnergy[0] +"\t\t"+ neighbourEnergy[1]);
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
		int searchStepsTotal = 3000;
		int searchStepsCicle = 4;

		double cambioPhi = 3.0;
		double cambioPsi = 3.0;

		String dirStartStruc = args[0];
		String dirTargetStruc = args[1];
		String fileDir = "out/";
		long initSeed = (long)(1000*Math.random());
		if(args.length>2){
			fileDir = args[2];
		}
		SimulatedAnnealing3DProtFromRCC simA = new SimulatedAnnealing3DProtFromRCC(searchStepsTotal,searchStepsCicle,cambioPhi,
																									cambioPsi,dirStartStruc,dirTargetStruc,fileDir,initSeed);
		//SimulatedAnnealing3DProtFromRCC simA = new SimulatedAnnealing3DProtFromRCC(args[0], args[1]);
		simA.initialize(2);
		simA.run(2);
	}
}
