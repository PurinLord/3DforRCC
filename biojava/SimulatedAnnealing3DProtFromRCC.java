import java.lang.Math;
import java.io.*;

import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.Chain;

import rccto3d.*;
import rccto3d.Trans;
import rccto3d.PDBfromFASTA;

public class SimulatedAnnealing3DProtFromRCC 
{

  // Set initial temp
	static double temp = 10000;

  // Cooling rate
  static double coolingRate = 0.003;

	//factor de cambio para Phi y Psi
	static double cambioPhi = 10.0;
	static double cambioPsi = 10.0;

 // Calculate the acceptance probability
 public static double acceptanceProbability(double energy, double newEnergy, double temperature) 
  {
   // If the new solution is better, accept it
   if (newEnergy < energy) return 1.0;
   // If the new solution is worse, calculate an acceptance probability
   return Math.exp((energy - newEnergy) / temperature);
  }

 public static int[] calcRCC(String s1){
	int rcc[] = new int[26];
	try{
	String s2 = "A";
	Process p = Runtime.getRuntime().exec("python create_26dvRCC.py "+s1+" "+s2);
	BufferedReader in = new BufferedReader(new InputStreamReader(p.getInputStream()));
	String ret = in.readLine();
	String val = "";
	int j = 0;
	for(int i = 0; i < ret.length(); i++){
		char c = ret.charAt(i);
		if(c==','){
			rcc[j] = Integer.valueOf(val);
			val = "";
			j++;
		}else{
			val+=c;
		}
	}
	rcc[j] = Integer.valueOf(val);
	}catch(Exception e){
		System.out.println("error" + e.toString());
	}
	return rcc;
}

 public static int[] calcRCC(Structure struc)
  {
   int[] result = new int[26];
   int i=0;

   for(i=0; i<26; i++) result[i]=0;

   return result;
  } // end calcRCC

	public static double calcSimilarity(int[] rcc1, int[] rcc2)
  {
		int sum = 0;
    int i=0;
		for(i=0; i<26; i++){
			sum += (rcc1[i] - rcc2[i])^2;
		}

   	return Math.sqrt(sum);
  }// end calcSimilarity

	public static void alterconformation(Structure struc){
		int largo = struc.size();
		//for(i=0; i<largo; i++){
		//}
		int dAmino = (int) (largo * Math.random());
		double dPsi = cambioPsi * Math.random();
		double dPhi = cambioPhi * Math.random();
		Chain chain = struc.getChain(0);
		Trans.rotatePsi(chain, dAmino, dPsi);
		Trans.rotatePhi(chain, dAmino, dPhi);
	}

 public static void main(String[] args) 
  {
   // Load Protein Sequence
		PDBfromFASTA pff = new PDBfromFASTA();
		String pdb = null;
		Structure struc_ini = null;

		try{
		pdb = pff.pdbFromFile(args[0], args[1]);
   //ProteinSequence seq = new ProteinsSequence(args[0]);

   // Build initial Protein 3D structure
		struc_ini = pff.shapeProtein(pdb);
		}catch(Exception e){
			System.out.println(e.toString());
		}
		Structure struc_model = null;
   //Protein3D prot3d_ini = new Protein3D(seq,args[1]);
   //Protein3D prot3d_model = null;
		PDBfromFASTA.writePDB("out/struc_ini.pdb", struc_ini);

   // Declaration of variables to store energy values
   double currentEnergy = 0, neighbourEnergy = 0;

   // Initialize intial solution
   double distance_ini = calcSimilarity(calcRCC(args[2]), calcRCC("out/struc_ini.pdb"));
   System.out.println("Initial solution distance: " + distance_ini);
   double distance_neighbor=0.0;

   // Set as current best
   double best = distance_ini;
        
   // Loop until system has cooled
   while (temp > 1) 
    {
     // Create new neighbour 3d model
			struc_model = PDBfromFASTA.readPDB("out/struc_ini.pdb");
     //seq_model = seq_ini.clone();

     // Get a random conformation for this new neighbor
     alterconformation(struc_model);
			PDBfromFASTA.writePDB("out/struc_model.pdb", struc_model);
            
     // Get energy of solutions
     currentEnergy = calcSimilarity(calcRCC(args[2]), calcRCC("out/struc_ini.pdb"));;
     neighbourEnergy = calcSimilarity(calcRCC(args[2]), calcRCC("out/struc_model.pdb"));;

     // Decide if we should accept the neighbour
     if (acceptanceProbability(currentEnergy, neighbourEnergy, temp) > Math.random()) {
       struc_ini = struc_model;
      }

     // Keep track of the best solution found
     if (neighbourEnergy < best) 
      {
       best = neighbourEnergy;
      }
            
     // Cool system
     temp *= 1-coolingRate;
    }

   System.out.println("Final solution distance: " + best);//.getDistance());
   System.out.println("Tour: " + best);
  }
}

