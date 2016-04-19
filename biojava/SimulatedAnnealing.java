package rccto3d.optimisation;

import java.util.Vector;
import java.util.Arrays;
import java.util.Random;
import java.lang.Math;
import java.io.*;

import rccto3d.*;

public class SimulatedAnnealing {
	// Set initial temp
	private double temp;
	private double tempPercent = 0.15;
	private int searchStepsTotal;
	private int searchStepsCicle;

	// Cooling rate
	private double coolingRate;

	private boolean rawAccept = false;

	private Random rdm;

	private StructureMannager sm;

	public SimulatedAnnealing(int searchStepsTotal, int searchStepsCicle, StructureMannager sm){
		this.searchStepsTotal = searchStepsTotal;
		this.searchStepsCicle = searchStepsCicle;
		this.temp = 0;
		this.coolingRate = 0;

		this.sm = sm;
		rdm =  new Random(System.currentTimeMillis());
	}

	public void setRawAccept(boolean rawAccept){
		this.rawAccept = rawAccept;
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
	
	public double calcInitialTemp(int verbos){

		int count = 0;
		int searchStepsCicle = 100;
		double sum = 0;
		double currentEnergy;
		while(count < searchStepsCicle){
			currentEnergy = sm.tempCalc();

			sum += currentEnergy;
			if(verbos > 2){
				System.out.println("energy " + currentEnergy + " sum " + sum);
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

	public void calcCoolRateStd(int verbos, int searchStepsTotal, int searchStepsCicle){
		double steps = searchStepsTotal/(double)searchStepsCicle;
		this.coolingRate = 1-(Math.pow((1/this.temp),(1.0/steps)));
	}
	
	public String initialize(int verbos){
		String outString = null;
		if(temp == 0){
			calcInitialTemp(verbos);
		}
		if(coolingRate == 0){
			calcCoolRateStd(verbos, searchStepsTotal, searchStepsCicle);
		}

		if (verbos > 0){
			outString = "temp "+temp +" coolRate "+coolingRate+" totalSearchSteps "+searchStepsTotal+" searchStepsCicle "+searchStepsCicle;
			if (verbos > 1){
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
		
		double currentEnergy[] = {0.0,0.0,0.0};
		double neighbourEnergy[] = {0.0,0.0,0.0};
		double distance_ini[] = {0.0,0.0,0.0};
		double best;

		// Initialize intial solution and set current as best
		distance_ini[0] = sm.calcEnergy();
		best = distance_ini[0];
		currentEnergy = Arrays.copyOf(distance_ini, distance_ini.length);
		
		if (verbos > 0){
			System.out.println("Initial solution distance: "+distance_ini[0]+" "+distance_ini[1]+" "+distance_ini[2]);
		}
		
		// Loop until system has cooled
		for (;temp > 1;temp = stdCooling(temp)){
			// Search steps
			sm.angleCooling();
			sm.addAngleStep(searchStepsCicle, searchStepsTotal);

			for(int step = 0; step < searchStepsCicle; step++){
				// Get a random conformation for this new neighbor
				sm.alterConformation();
				// Get energy of solution
				neighbourEnergy[0] = sm.calcEnergy();
			
				// Decide if we should accept the neighbour
				if (acceptanceProbability(currentEnergy[0], neighbourEnergy[0], temp) > rdm.nextDouble()) {
					sm.acceptModel();
					//PDBfromFASTA.writePDB(fileDir + "struc_fit.pdb", struc_fit);
					currentEnergy = Arrays.copyOf(neighbourEnergy, neighbourEnergy.length);
					if (verbos > 0){
						if (verbos > 1){
							long time = System.nanoTime() - elapsedTime;
							System.out.println(temp + "> " + currentEnergy[0] + " time: " + time/3600000000000.0);
						}else{
							System.out.println(temp + "> " + currentEnergy[0]);
						}
					}
				}
			
				// Keep track of the best solution found
				if (neighbourEnergy[0] < best){
					best = neighbourEnergy[0];
					sm.setFitAsBeast();
					//PDBfromFASTA.writePDB(fileDir + "sol_" + best + ".pdb", struc_fit);
					if(best == 0){
						temp = 0;
					}
				}
				if (verbos > 1){
					sm.printEnergy();
				}
			}
		}
		
		if (verbos > 0){
			System.out.println("Final solution distance: " + best);//.getDistance());
			elapsedTime = System.nanoTime() - elapsedTime;
			System.out.println("Total execution time: " + elapsedTime/3600000000000.0);
		}
		sm.saveBest();
		return best;
	}

	public static void main(String[] args){
		int searchStepsTotal = 3000;
		int searchStepsCicle = 1;
	
		double anguloInicial = 0.05;
		double anguloFinal = 0.001;
	
		String[] fileName = new String[4];
		fileName [0] = args[0];
		fileName [1] = args[1];
		fileName [2] = args[2];
		fileName [3] = args[3];
		String fileDir = "out/";
		long initSeed = (long)(1000*Math.random());
		//if(args.length>2){
		//	fileDir = args[2];
		//}

		StructureMannager sm = new StructureMannager(0, 1, 1, 2);
		sm.setConditions(anguloInicial, anguloFinal, fileDir, initSeed);
		sm.loadStructures(fileName, 2);

		SimulatedAnnealing siA = new SimulatedAnnealing(searchStepsTotal,searchStepsCicle,sm);
		siA.setTemp(2.5);
		siA.initialize(2);
		siA.run(2);
	}
}
