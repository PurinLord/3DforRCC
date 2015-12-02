package rccto3d.optimisation;

import java.io.File;

public class Refiner{

	static int searchStepsTotal;
	static int searchStepsCicle;
	static double anguloInicial;
	static double anguloFinal;
	static String dirStartStruc;
	static String dirTargetStruc;
	static String fileDir;
	static long initSeed;

	public Refiner(String dirStartStruc, String dirTargetStruc, String fileDir){
		searchStepsTotal = 1000;
		searchStepsCicle = 1;
		anguloInicial = 0.01;
		anguloFinal = 0.001;
		this.dirStartStruc = dirStartStruc;
		this.dirTargetStruc = dirTargetStruc;
		this.fileDir = fileDir;
		initSeed = (long)(1000*Math.random());
	}

	public static String simulA(){
		SimulatedAnnealing3DProtFromRCC simA = new SimulatedAnnealing3DProtFromRCC(searchStepsTotal,searchStepsCicle,anguloInicial,
																									anguloFinal,dirStartStruc,dirTargetStruc,fileDir,initSeed);
		simA.setTemp(2.5);
		simA.initialize(1);
		double best =  simA.run(0);
		return ("sol_" + best + ".pdb");
	}

	public static String simulA(String start){
		String dirStartStruc = "";
		try{
			File file = new File(start);
			if (file.exists()) {
				dirStartStruc = start;
			}else{
				dirStartStruc = fileDir+"struc_ini.pdb";
			}
		}catch(Exception e){
			e.printStackTrace();
		}
		SimulatedAnnealing3DProtFromRCC simA = new SimulatedAnnealing3DProtFromRCC(searchStepsTotal,searchStepsCicle,anguloInicial,
																									anguloFinal,dirStartStruc,dirTargetStruc,fileDir,initSeed);
		simA.setTemp(2.5);
		simA.initialize(1);
		double best =  simA.run(0);
		return ("sol_" + best + ".pdb");
	}

	public static void main(String[] args){
		Refiner ref = new Refiner(args[0], args[1], args[2]);
		String sol = ref.simulA();
		for(int i=1;i<20;i++){
			sol = ref.simulA(Refiner.fileDir+sol);
		}
	}
}
