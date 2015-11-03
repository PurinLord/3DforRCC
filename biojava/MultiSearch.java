package rccto3d.optimisation;

import java.lang.Integer;
import java.io.File;
import java.io.BufferedWriter;
import java.io.FileWriter;
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
	}

	public void run(){
		System.out.println(fileDir +" "+ simA.run(0));	
	}

	public void write(String filename, String string){
		try{
			File file = new File(fileDir+filename);
			if (!file.exists()) {
				file.createNewFile();
			}
			synchronized(file) {
				BufferedWriter bw = new BufferedWriter(new FileWriter(file.getAbsoluteFile()));
				bw.write(string);
				bw.close();
			}
		}catch(Exception e){
			e.printStackTrace();
		}
	}

	public static void main(String args[]){

		int searchStepsTotal = 3500;
		int searchSteps = 2;
		double cambioPhi = 36.0;
		double cambioPsi = 36.0;
		String dirStartStruc = args[0];
		String dirTargetStruc = args[1];
		String fileDir = "out/";
		double temp = 0;

		MultiSearch mSearch;
		mSearch = new MultiSearch(searchStepsTotal,searchSteps,cambioPhi,cambioPsi,dirStartStruc,dirTargetStruc,
																fileDir + 0 + "/",0);
		mSearch.simA.initialize(0);
		temp = mSearch.simA.getTemp();
		for(int i = 1; i < Integer.parseInt(args[2]); i++){
			searchStepsTotal = 3500;
			searchSteps = 2;
			cambioPhi = 36.0;
			cambioPsi = 36.0;

			mSearch = new MultiSearch(searchStepsTotal,searchSteps,cambioPhi,cambioPsi,dirStartStruc,dirTargetStruc,
																	fileDir + i + "/",i);
			mSearch.simA.setTemp(temp);
			mSearch.write("rcc.our", mSearch.simA.initialize(1));

			(new Thread(mSearch)).start();
		}
	}
}
