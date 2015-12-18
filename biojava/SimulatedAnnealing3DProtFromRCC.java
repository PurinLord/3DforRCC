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
	private double initPhi;
	private double initPsi;
	private double minPhi;
	private double minPsi;
	private double cambioPhi;
	private double cambioPsi;
	private double angleSteps = 0;

	private String dirStartStruc;
	private String fastaID = "a";
	private String dirTargetStruc;
	private String fileDir = "";

	//private int divNumSegment;
	//private int divMaxSize;
	//private int divMinSize;
	//private int divUndefMax;
	//private int divUndefMin;

	private long initSeed;
	// 0 - RCC
	// 1 - RMSD
	private int energyType = 1;
	private boolean dualEnergy = false;
	private boolean rawAccept = false;

	private Substitutor sub = null;
	private Structure struc_subSeed = null;

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

		initPhi = cambioPhi = 9.0;
		initPsi = cambioPsi = 9.0;
		this.minPhi = 0.01;
		this.minPsi = 0.01;
		this.dirStartStruc = dirStartStruc;
		this.dirTargetStruc = dirTargetStruc;
		fileDir = "out/";
		initSeed = (long)(1000*Math.random());
	}

	public SimulatedAnnealing3DProtFromRCC(int searchStepsTotal, int searchStepsCicle, double anguloInicial, double anguloFinal,
																					String dirStartStruc, String dirTargetStruc, String fileDir, long initSeed){
		this.searchStepsTotal = searchStepsTotal;
		this.searchStepsCicle = searchStepsCicle;
		this.temp = 0;
		this.coolingRate = 0;

		this.initPhi = this.cambioPhi = this.initPsi = this.cambioPsi = anguloInicial;
		this.minPhi = this.minPsi = anguloFinal;

		this.dirStartStruc = dirStartStruc;
		this.dirTargetStruc = dirTargetStruc;
		this.fileDir = fileDir;
		this.initSeed = initSeed;
	}

	public void setEnergyType(int type){
		this.energyType = type;
	}

	public void setRawAccept(boolean rawAccept){
		this.rawAccept = rawAccept;
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

	public void setSubsitutor(Substitutor sub, Structure struc_subSeed){
		this.sub = sub;
		this.struc_subSeed = struc_subSeed;
	}

	public void setSubsitutor(Structure struc_subSeed){
		this.struc_subSeed = struc_subSeed;
	}

	public void setSubsitutor(String dirStrucSeed){
		Structure struc_seed =PDBfromFASTA.readPDB(dirStrucSeed);
		this.struc_subSeed = struc_subSeed;
	}

	//Calculate the acceptance probability
	public double acceptanceProbability(double energy, double newEnergy, double temperature) {
		if(rawAccept){
			return rawAcceptanceProbability(energy, newEnergy);
   	}
		// If the new solution is better, accept it
		if (newEnergy < energy) return 1.0;
		// If the new solution is worse, calculate an acceptance probability
		return Math.exp((energy - newEnergy) / temperature);
	}

	public double rawAcceptanceProbability(double energy, double newEnergy) {
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

	public double calcRMSD3(Structure struc_current, Structure struc_target){
		double currentRMSd = 0.0;
		try {
 			StructureAlignment algorithm  = StructureAlignmentFactory.getAlgorithm(CeMain.algorithmName);
 			
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
			energy[1] = calcRMSD3(s1, s2);
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
			energy[1] = calcRMSD3(s1, s2);
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

	public Structure alterConformationParts(Structure struc){
		Structure alterStruc = (Structure)struc.clone();
		Chain chain = alterStruc.getChain(0);
		Random rdm =  new Random(System.currentTimeMillis());
		int largo = chain.getAtomLength();
		double dPsi;
		double dPhi;
		int index = 0;
		Vector<Integer> segment = sub.getDivition().elementAt(index);
		for(int i=0; i<largo; i++){
			if(i == segment.elementAt(0)){
				i += segment.elementAt(1);
				index++;
				if(index < sub.getDivition().size()){
					segment = sub.getDivition().elementAt(index);
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

	//
	//Random de 0 error
	//al dividir salen 0 donde no deberían
	//todo mayor que 0!!!!!
	public Substitutor createSubstitutor(Structure subStructure){
		Substitutor sub = new Substitutor(subStructure);
		Random rdm =  new Random(System.currentTimeMillis());
		int largo = subStructure.getChain(0).getAtomLength();
		int numSegment, maxSize, minSize, undefMax, undefMin;
		numSegment = rdm.nextInt(largo) + 2;
		maxSize = rdm.nextInt(largo) + 5;
		undefMax = rdm.nextInt(largo) + 1;
		if(undefMax > maxSize){int tmp = maxSize; maxSize = undefMax; undefMax = tmp;}//SWAP
		minSize = rdm.nextInt(maxSize);
		undefMin = rdm.nextInt(undefMax);
		int largoActual = numSegment*maxSize + numSegment*undefMax + undefMin;
		numSegment *= Math.sqrt(largo/(float)largoActual);
		maxSize *= Math.sqrt(largo/(float)largoActual);
		undefMax *= Math.sqrt(largo/(float)largoActual);
		undefMin *= Math.sqrt(largo/(float)largoActual);
		minSize *= Math.sqrt(largo/(float)largoActual);
		if(maxSize == 0){maxSize = 1;}
		if(undefMax == 0){undefMax = 1;}
		if(minSize == 0){minSize = 1;}
		if(undefMin == 0){undefMin = 1;}
		largoActual = numSegment*maxSize + numSegment*undefMax + undefMin;
		//System.out.println(numSegment+"\t"+maxSize+"\t"+minSize+"\t"+undefMax+"\t"+undefMin+"\t"+largoActual);

		sub.createDivition(numSegment, maxSize, minSize, undefMax, undefMin);
		return sub;
	}

	public double angleCooling(double time){
		return initPhi*Math.exp(-time/(-1.0/Math.log(minPhi/initPhi)));
	}

	public double stdCooling(double temp){
		return temp * (1-coolingRate);
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
	//			 3 - temp calc
	
	public double calcInitialTemp(int verbos, Structure struc_subSeed){
		PDBfromFASTA.writePDB(fileDir + "calcTmp0.pdb", struc_subSeed);
		int seedRCC[] = calcRCC(fileDir + "calcTmp0.pdb");
		//int alterRCC[];
		Structure alterStruc = (Structure)struc_subSeed.clone();
		int count = 0;
		int searchStepsCicle = 100;
		//double max = 0;
		double sum = 0;
		double currentEnergy[]={0.0,0.0};
		while(count < searchStepsCicle){
			alterStruc = alterConformationAll(alterStruc, 180, 180);

			currentEnergy = calcEnergy(struc_subSeed,alterStruc, seedRCC);
			//PDBfromFASTA.writePDB(fileDir + "calcTmp1.pdb", alterStruc);
			//alterRCC = calcRCC(fileDir + "calcTmp1.pdb");
			//currentEnergy = calcSimilarity(seedRCC, alterRCC);;
			//currentEnergy[1] = calcRMSD3(struc_subSeed,alterStruc);

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
		if(energyType == 0 || dualEnergy){
			targetRCC = calcRCC(dirTargetStruc);
		}

		if(struc_subSeed != null){
			if(sub == null){
				sub = createSubstitutor(struc_ini);
			}else{
				struc_ini = sub.fakeSubstitute(struc_subSeed);
			}
		}

		if(temp == 0){
			calcInitialTemp(verbos, struc_target);
		}
		if(coolingRate == 0){
			calcCoolRateStd(verbos, searchStepsTotal, searchStepsCicle);
		}

		if (verbos > 0){
			outString = dirStartStruc + " " + dirTargetStruc + 
				"\ntemp "+temp +" coolRate "+coolingRate+" totalSearchSteps "+searchStepsTotal+" searchStepsCicle "+searchStepsCicle +
				"\nangulo Inicial " + initPhi + " angulo final " + minPhi + " initSeed " + initSeed;
			if(struc_subSeed != null){
				Vector<Vector<Integer>> divition = sub.getDivition();
				outString += "\n" + divition;
			}
			if (verbos > 1){
				if(energyType == 0 || dualEnergy){
					for(int i : targetRCC){System.out.print(i + " ");}
				}
				System.out.println(outString);
				System.out.print("\n");
			}
			return outString;
		}
		return null;
	}

	public double run(int verbos){
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
		if(energyType == 0 || dualEnergy){
			currentRCC = calcRCC(fileDir + "struc_ini.pdb");
			distance_ini[0] = calcSimilarity(targetRCC, currentRCC);
		}
		distance_ini[1] = calcRMSD3(struc_ini,struc_target);
		
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
			cambioPhi = angleCooling(angleSteps);
			cambioPsi = cambioPhi;
			angleSteps += (double)searchStepsCicle/searchStepsTotal;
			for(int step = 0; step < searchStepsCicle; step++){
				// Get a random conformation for this new neighbor
				if(struc_subSeed != null){
					struc_model = alterConformationParts(struc_fit);
				}else{
					if(target != null){
						struc_model = alterConformationAll(struc_fit, target);
					}else{
						struc_model = alterConformationAll(struc_fit);
					}
				}
				// Get energy of solution
				neighbourEnergy = calcEnergy(struc_model,struc_target , targetRCC);
				//PDBfromFASTA.writePDB(fileDir + "struc_model.pdb", struc_model);
				//modelRCC = calcRCC(fileDir + "struc_2.pdb");
				//neighbourEnergy[0] = calcSimilarity(targetRCC, modelRCC);
				//neighbourEnergy[1] = calcRMSD3(struc_target, struc_model);
			
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
					//if(energyType == 0 || dualEnergy){
					//	for(int i : modelRCC){System.out.print(i + " ");}
					//	System.out.print("\n");
					//}
					if(dualEnergy){
						System.out.println(neighbourEnergy[0] +"\t\t"+ neighbourEnergy[1]);
					}else{
						System.out.println(neighbourEnergy[energyType] +"\t\t"+ cambioPhi);
					}
				}
			}
		}
		
		if (verbos > 0){
			System.out.println("Final solution distance: " + best);//.getDistance());
			elapsedTime = System.nanoTime() - elapsedTime;
			System.out.println("Total execution time: " + elapsedTime/3600000000000.0);
		}
		return best;
  }

	public static void main(String[] args){
		int searchStepsTotal = 3000;
		int searchStepsCicle = 1;

		double anguloInicial = 0.05;
		double anguloFinal = 0.001;

		String dirStartStruc = args[0];
		String dirTargetStruc = args[1];
		String fileDir = "out/";
		long initSeed = (long)(1000*Math.random());
		if(args.length>2){
			fileDir = args[2];
		}
		SimulatedAnnealing3DProtFromRCC simA = new SimulatedAnnealing3DProtFromRCC(searchStepsTotal,searchStepsCicle,anguloInicial,
																									anguloFinal,dirStartStruc,dirTargetStruc,fileDir,initSeed);
		//SimulatedAnnealing3DProtFromRCC simA = new SimulatedAnnealing3DProtFromRCC(args[0], args[1]);
		simA.setTemp(2.5);
		//Substitutor sub = new Substitutor(Trans.readPDB(args[1]));
		//sub.createDivition(4, 20, 20, 4, 2);
		simA.setSubsitutor(PDBfromFASTA.readPDB(args[3]));
		simA.initialize(2);
		simA.run(2);
	}
}
