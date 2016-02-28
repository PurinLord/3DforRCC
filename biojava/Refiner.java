package rccto3d.optimisation;

import java.util.*;
import java.io.File;
import java.io.BufferedReader;
import java.io.FileReader;
import java.util.Vector;
import java.lang.Double;
import java.util.PriorityQueue;

import org.biojava.bio.structure.StructureImpl;

import rccto3d.*;

public class Refiner implements Runnable{
	int searchStepsTotal;
	int searchStepsCicle;
	double anguloInicial;
	double anguloFinal;
	double iniFin;
	double finalFin;
	double iniPend;
	double finalPend;
	int numBusquedas;
	String startDir;
	String dirStartStruc;
	String dirTargetStruc;
	String fileDir;
	long initSeed;
	double beast;

	public Refiner(String rootDir, int startIndex, double startValue, String targetStruc, int i){
		searchStepsTotal = 2000;
		searchStepsCicle = 1;
		anguloInicial = 0.01;
		anguloFinal = 0.001;
		iniFin = 0.0001;
		finalFin = 0.000001;
		numBusquedas = 40;
		iniPend = (anguloInicial-iniFin)/numBusquedas;
		finalPend = (anguloFinal-finalFin)/numBusquedas;
		this.startDir = rootDir+startIndex+"/";
		this.dirStartStruc = startDir+"sol_"+startValue+".pdb";
		this.dirTargetStruc = targetStruc;
		this.fileDir = rootDir+"ref_"+(i+1)+"/";
		initSeed = (long)(1000*Math.random());
		beast = 0;
	}

	public Substitutor setSubsitutor(String startDir){
 		String lastLine = "";
		Substitutor sub = new Substitutor(new StructureImpl());
		try{
			File file = new File(startDir+"rep.out");
			BufferedReader br = new BufferedReader(new FileReader(file));
			String line = "";
    	while ((line = br.readLine()) != null){ lastLine=line;}
			sub.readDivition(lastLine);
		}catch(Exception e){
			e.printStackTrace();
		}
		return sub;
	}

	public void run(){
		File f = new File(fileDir);
		if (!f.exists()){
			f.mkdirs();
		}
		for(int i=numBusquedas;i>0;i--){
			try{
				File file = new File(dirStartStruc);
				if (!file.exists()) {
					dirStartStruc = fileDir+"struc_ini.pdb";
			}}catch(Exception e){
				e.printStackTrace();
			}
			//System.out.println(searchStepsTotal+" "+searchStepsCicle+" "+anguloInicial+" "+anguloFinal+" "+dirStartStruc+" "+dirTargetStruc+" "+fileDir+" "+initSeed);
			SimulatedAnnealing3DProtFromRCC simA = new SimulatedAnnealing3DProtFromRCC(searchStepsTotal,
					searchStepsCicle,anguloInicial,anguloFinal,dirStartStruc,dirTargetStruc,fileDir,initSeed);
			simA.setRawAccept(true);
			simA.setEnergyType(1);
			simA.setTemp(4);
			simA.initialize(0);
			simA.setSubsitutor(this.setSubsitutor(startDir));
			beast = simA.run(0);
			dirStartStruc = fileDir + "sol_" + beast + ".pdb";
			//System.out.println(beast);
			anguloInicial = iniFin+(i*iniPend);
			anguloFinal = finalFin+(i*finalPend);
		}
		System.out.println(fileDir +" "+ beast);	
	}

	public static void main(String[] args){
		String rootDir = args[0];
		String reportFileName = "rep.out";
		String targetStruc = rootDir+"1/struc_target.pdb";
		ScanReport scan = new ScanReport();
		Vector<Pair> vec = scan.scan(rootDir+reportFileName, Integer.parseInt(args[1]));
		Refiner ref;
		for(int i=0;i<vec.size();i++){
			double startValue = vec.elementAt(i).energy;
			int startIndex = vec.elementAt(i).searchNum;
			ref = new Refiner(rootDir,startIndex,startValue, targetStruc,i);
			//ref = new Refiner(rootDir+startDir+"/sol_"+startName+".pdb", targetStruc, rootDir+"ref_"+(i+1)+"/");
			//ref.refine();
			(new Thread(ref)).start();
		}
	}
}

class ScanReport{
	public Vector<Pair> scan (String referencia, int numRefine){
		Vector<Pair> refineGrup = new Vector<Pair>(numRefine);
		try {
			BufferedReader br = new BufferedReader(new FileReader(referencia)); 
			String line;
			int searchNum;
			double energy;
			PriorityQueue<Pair> queue = new PriorityQueue<Pair>();
			Pair pair;
			while ((line = br.readLine()) != null) {
				energy = Double.parseDouble(line.substring(line.lastIndexOf("/") + 1));
				searchNum = Integer.parseInt(line.substring(line.indexOf("/") + 1, line.lastIndexOf("/")));
				pair = new Pair(searchNum, energy);
				queue.add(pair);
				//System.out.println(line +" "+ pair.energy +" "+ pair.searchNum);
			}
			for(int i=0;i<numRefine;i++){
				if((pair = queue.poll()) != null){
					refineGrup.add(pair);
				}else{
					break;
				}
			}
		}catch(Exception e){
			e.printStackTrace();
		}
		return refineGrup;
	}
}

class Pair implements Comparable<Pair>{
	int searchNum;
	double energy;
	
	Pair(int searchNum, double energy){
		this.searchNum = searchNum;
		this.energy = energy;
	}

	public int compareTo(Pair p){
		if(this.energy > p.energy) return 1;
		return -1;
	}
}
