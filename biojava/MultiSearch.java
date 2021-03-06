package rccto3d.optimisation;

import java.lang.Integer;
import java.util.Random;
import java.io.File;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.lang.Thread;

public class MultiSearch implements Runnable{

	private SimulatedAnnealing siA;
	private StructureMannager sm;
	private String fileDir;
	String reporte;

	public MultiSearch(int searchStepsTotal, int searchSteps, double angIni, double angFin,
			String[] fileName, String fileDir, long initSeed){

		this.fileDir = fileDir;
		File f = new File(fileDir);
		if (!f.exists()){
			f.mkdirs();
		}

		sm = new StructureMannager(0, 0, 0, 3);
		sm.setConditions(angIni, angFin, fileDir, initSeed);
		this.reporte = sm.loadStructures(fileName, 1);
		
		this.siA = new SimulatedAnnealing(searchStepsTotal, searchSteps, sm);
	}

	public void run(){
		System.out.println(fileDir +" "+ siA.run(0));	
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
		int searchSteps = 5;
		double anguloInicial = 9;
		double anguloFinal = 0.001;
		String[] fileName = new String[args.length - 1];
		fileName[0] = args[0];
		fileName[1] = args[1];
		String fileDir = "out/";
		double temp = 0;

		int threadSize = 4;
		int defaultThreads = Thread.getAllStackTraces().size();
		int liveThreads = 0;

		MultiSearch mSearch;
		//mSearch.simA.initialize(0);
		//temp = mSearch.simA.getTemp();
		temp = 2.5;
		Random rdm = null;
		int i = 1;
		for(i = 1; i <= threadSize; i++){
			rdm =  new Random(System.currentTimeMillis());
			mSearch = new MultiSearch(searchStepsTotal,searchSteps,anguloInicial,anguloFinal,fileName,fileDir + i + "/",i*rdm.nextLong());
			mSearch.siA.setTemp(temp);
			mSearch.write("rep.out", mSearch.siA.initialize(1)+"\n"+mSearch.reporte);
			(new Thread(mSearch)).start();
		}
		liveThreads = Thread.getAllStackTraces().size();
		for(; i <= Integer.parseInt(args[2]); i++){
			while(liveThreads - defaultThreads >= threadSize){
				liveThreads = Thread.getAllStackTraces().size();
			}

			rdm =  new Random(System.currentTimeMillis());
			mSearch = new MultiSearch(searchStepsTotal,searchSteps,anguloInicial,anguloFinal,fileName,fileDir + i + "/",i*rdm.nextLong());
			mSearch.siA.setTemp(temp);
			mSearch.write("rep.out", mSearch.siA.initialize(1)+"\n"+mSearch.reporte);
			(new Thread(mSearch)).start();

			liveThreads = Thread.getAllStackTraces().size();
		}
	}
}
