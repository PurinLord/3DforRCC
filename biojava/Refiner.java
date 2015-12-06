package rccto3d.optimisation;

import java.io.File;

public class Refiner{

	static int searchStepsTotal;
	static int searchStepsCicle;
	static double anguloInicial;
	static double anguloFinal;
	static double iniFin;
	static double finalFin;
	static double iniPend;
	static double finalPend;
	static int numBusquedas;
	static String dirStartStruc;
	static String dirTargetStruc;
	static String fileDir;
	static long initSeed;
	static double beast;

	public Refiner(String dirStartStruc, String dirTargetStruc, String fileDir){
		searchStepsTotal = 20000;
		searchStepsCicle = 1;
		anguloInicial = 0.01;
		anguloFinal = 0.001;
		iniFin = 0.0001;
		finalFin = 0.000001;
		numBusquedas = 40;
		iniPend = (anguloInicial-iniFin)/numBusquedas;
		finalPend = (anguloFinal-finalFin)/numBusquedas;
		this.dirStartStruc = dirStartStruc;
		this.dirTargetStruc = dirTargetStruc;
		this.fileDir = fileDir;
		initSeed = (long)(1000*Math.random());
		beast = 0;
	}

	public static double simulA(){
		try{
			File file = new File(dirStartStruc);
			if (!file.exists()) {
				dirStartStruc = fileDir+"struc_ini.pdb";
			}
		}catch(Exception e){
			e.printStackTrace();
		}
		SimulatedAnnealing3DProtFromRCC simA = new SimulatedAnnealing3DProtFromRCC(searchStepsTotal,searchStepsCicle,anguloInicial,
																									anguloFinal,dirStartStruc,dirTargetStruc,fileDir,initSeed);
		simA.setEnergyType(1);
		simA.setTemp(4);
		simA.initialize(1);
		beast = simA.run(0);
		return beast;
	}

	public static void main(String[] args){
		Refiner ref = new Refiner(args[0], args[1], args[2]);
		for(int i=ref.numBusquedas;i>0;i--){
			ref.dirStartStruc = ref.fileDir + "sol_" + ref.simulA() + ".pdb";
			System.out.println(ref.beast);
			ref.anguloInicial =  ref.iniFin+(i*ref.iniPend);
			ref.anguloFinal = ref.finalFin+(i*ref.finalPend);
		}
	}
}
