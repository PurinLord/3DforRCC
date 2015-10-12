package rccto3d.optimisation;

import rccto3d.PDBfromFASTA;

import java.lang.Integer;
import java.io.File;
import java.lang.Thread;

public class MultiSearch implements Runnable{

	private SimulatedAnnealing3DProtFromRCC simA;
	private String fileDir;

	public MultiSearch(double temp, int searchSteps, double coolingRate, double cambioPhi, double cambioPsi, String startPdb,
											int targetRCC[], String fileDir, long initSeed){
		
		this.simA = new SimulatedAnnealing3DProtFromRCC(temp, searchSteps, coolingRate, cambioPhi, 
																																	cambioPsi, startPdb, targetRCC, fileDir, initSeed);
		//SimulatedAnnealing3DProtFromRCC simA = new SimulatedAnnealing3DProtFromRCC(dirStartStruc, dirTargetStruc);
		this.fileDir = fileDir;
	}

	public void run(){
		File f = new File(fileDir);
		if (!f.exists()){
			boolean success = f.mkdirs();
		}
		simA.run(0);	
	}

	public static void main(String args[]){

		double temp = 10000;
		int searchSteps = 2;
  	double coolingRate = 0.003;
		double cambioPhi = 36.0;
		double cambioPsi = 36.0;
		String dirStartStruc = args[0];
		String dirTargetStruc = args[1];
		String fileDir = "out/";
		long initSeed = 9;

		MultiSearch mSearch;
		//TODO (make precreation of rcc and pdb reazonable)
		PDBfromFASTA pff = new PDBfromFASTA();
		String startPdb = pff.pdbFromFile(dirStartStruc, "a");
		int targetRCC[] = SimulatedAnnealing3DProtFromRCC.calcRCC(dirTargetStruc, fileDir);
		for(int i = 1; i <= Integer.parseInt(args[2]); i++){
			mSearch = new MultiSearch(temp,searchSteps,coolingRate,cambioPhi,cambioPsi,startPdb,targetRCC,
																	fileDir + i + "/",i);
			(new Thread(mSearch)).start();
		}
	}
}
