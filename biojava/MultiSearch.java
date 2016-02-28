package rccto3d.optimisation;

import java.lang.Integer;
import java.util.Random;
import java.io.File;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.lang.Thread;
import java.util.concurrent.CyclicBarrier;

public class MultiSearch implements Runnable{

	protected SimulatedAnnealing3DProtFromRCC simA;
	protected String fileDir;
	CyclicBarrier cb = null;

	public MultiSearch(int searchStepsTotal, int searchSteps, double angIni, double angFin,
			String dirStartStruc,String dirTargetStruc, String fileDir, long initSeed, CyclicBarrier cb){
		
		this.cb = cb;
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
		try{
			double result = simA.run(0);
			cb.await();
			System.out.println(fileDir +" "+ result);	
		}catch(Exception e){
			e.printStackTrace();
		}
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
		Group g = new Group(Integer.parseInt(args[2]), 8);
		g.run(args[0], args[1]);
	}
}

class Group{
	int total = 0;
	int size = 0;
	CyclicBarrier cbl = null;
	public Group(int total, int size){
		this.total = total/size;
		this.size = size;
		cbl = new CyclicBarrier(size);
	}

	public void run(String dirStart, String dirTarget){
		Random rdm =  new Random(System.currentTimeMillis());
		int searchStepsTotal = 3000;
		int searchSteps = 3;
		double anguloInicial = 9;
		double anguloFinal = 0.001;
		String dirStartStruc = dirStart;
		String dirTargetStruc = dirTarget;
		String fileDir = "out/";
		double temp = 0;

		MultiSearch mSearch;
		temp = 2.5;
		for(int i = 0; i < total; i++){
			for(int j = 1; j <= size; j++){

				mSearch = new MultiSearch(searchStepsTotal,searchSteps,anguloInicial,anguloFinal,dirStartStruc,
						dirTargetStruc,fileDir + ((i*size)+j) + "/",i*rdm.nextLong(), cbl);
				mSearch.simA.setTemp(temp);
				mSearch.simA.setSubsitutor(dirTarget);
				mSearch.write("rep.out", mSearch.simA.initialize(1));

				(new Thread(mSearch)).start();
			}
		}
	}
}
