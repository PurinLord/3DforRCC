package rccto3d.optimisation;

import java.lang.Integer;
import java.io.File;
import java.lang.Thread;

public class MultiSearch implements Runnable{

	private SimulatedAnnealing3DProtFromRCC simA;
	private String fileDir;

	public MultiSearch(int searchStepsTotal, int searchSteps, double cambioPhi, double cambioPsi, String dirStartStruc,
											String dirTargetStruc, String fileDir, long initSeed){
		
		this.simA = new SimulatedAnnealing3DProtFromRCC(searchStepsTotal, searchSteps, cambioPhi, cambioPsi, dirStartStruc,
																																	 dirTargetStruc, fileDir, initSeed);
		//SimulatedAnnealing3DProtFromRCC simA = new SimulatedAnnealing3DProtFromRCC(dirStartStruc, dirTargetStruc);
		this.fileDir = fileDir;
		File f = new File(fileDir);
		if (!f.exists()){
			f.mkdirs();
		}
		this.simA.initialize(0);
	}

	public void run(){
		System.out.println(fileDir +" "+ simA.run(0));	
	}

	public static void main(String args[]){

		int searchStepsTotal = 3500;
		int searchSteps = 2;
		double cambioPhi = 36.0;
		double cambioPsi = 36.0;
		String dirStartStruc = args[0];
		String dirTargetStruc = args[1];
		String fileDir = "out/";
		long initSeed = 9;

		MultiSearch mSearch;
		for(int i = 1; i <= Integer.parseInt(args[2]); i++){
			mSearch = new MultiSearch(searchStepsTotal,searchSteps,cambioPhi,cambioPsi,dirStartStruc,dirTargetStruc,
																	fileDir + i + "/",i);

			(new Thread(mSearch)).start();
		}
	}
}
