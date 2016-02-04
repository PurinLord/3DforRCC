package rccto3d.optimisation;

import java.lang.Integer;
import java.util.Random;
import java.io.File;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.lang.Thread;

public class MultiSearch implements Runnable{

	private SimulatedAnnealing3DProtFromRCC simA;
	private String fileDir;

	public MultiSearch(int searchStepsTotal, int searchSteps, double angIni, double angFin,
			String dirStartStruc,String dirTargetStruc, String fileDir, long initSeed){
		
		this.simA = new SimulatedAnnealing3DProtFromRCC(searchStepsTotal, searchSteps, angIni, angFin,
				dirStartStruc,dirTargetStruc, fileDir, initSeed);
		//SimulatedAnnealing3DProtFromRCC simA = new SimulatedAnnealing3DProtFromRCC(dirStartStruc, dirTargetStruc);
		this.simA.setEnergyType(1);
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
		int searchStepsTotal = 3000;
		int searchSteps = 3;
		double anguloInicial = 9;
		double anguloFinal = 0.001;
		String dirStartStruc = args[0];
		String dirTargetStruc = args[1];
		String fileDir = "out/";
		double temp = 0;

		MultiSearch mSearch;
		//mSearch = new MultiSearch(searchStepsTotal,searchSteps,cambioPhi,cambioPsi,dirStartStruc,dirTargetStruc,
		//														fileDir + 0 + "/",0);
		//mSearch.simA.initialize(0);
		//temp = mSearch.simA.getTemp();
		temp = 2.5;
		Random rdm = null;
		long elapsedTime = 0;
		elapsedTime = System.nanoTime();
		for(int i = 1; i <= Integer.parseInt(args[2]); i++){
			rdm =  new Random(System.currentTimeMillis());

			mSearch = new MultiSearch(searchStepsTotal,searchSteps,anguloInicial,anguloFinal,dirStartStruc,
					dirTargetStruc,fileDir + i + "/",i*rdm.nextLong());
			mSearch.simA.setTemp(temp);
			mSearch.simA.setSubsitutor(args[1]);
			mSearch.write("rep.out", mSearch.simA.initialize(1));

			(new Thread(mSearch)).start();
		}
		long time = System.nanoTime() - elapsedTime;
		System.out.println(temp + "> time: " + time/3600000000000.0);
	}
}
