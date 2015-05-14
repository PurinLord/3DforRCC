public class SimulatedAnnealing3DProtFromRCC 
{
 // Calculate the acceptance probability
 public static double acceptanceProbability(int energy, int newEnergy, double temperature) 
  {
   // If the new solution is better, accept it
   if (newEnergy < energy) return 1.0;
   // If the new solution is worse, calculate an acceptance probability
   return Math.exp((energy - newEnergy) / temperature);
  }

 public static int[] calcRCC(String f)
  {
   int[] result = new int[26];
   int i=0;

   for(i=0; i<26; i++) result[i]=0;

   return result;
  } // end calcRCC

 public static int[] calcRCC(Protein3D prot3d)
  {
   int[] result = new int[26];
   int i=0;

   for(i=0; i<26; i++) result[i]=0;

   return result;
  } // end calcRCC

 public static double calcSimilarity(int[] rcc1, int[] rcc2)
  {
   double result=0.0;


   return result;
  }// end calcSimilarity

 public static void main(String[] args) 
  {
   // Load Protein Sequence
   ProteinSequence seq = new ProteinsSequence(args[0]);

   // Build initial Protein 3D structure
   Protein3D prot3d_ini = new Protein3D(seq,args[1]);
   Protein3D prot3d_model = null;

   // Set initial temp
   double temp = 10000;

   // Cooling rate
   double coolingRate = 0.003;

   // Declaration of variables to store energy values
   int currentEnergy = 0, neighbourEnergy = 0;

   // Initialize intial solution
   double distance_ini = calcSimilarity(calcRCC(args[2]), calcRCC(prot3d_ini));
   System.out.println("Initial solution distance: " + distance_ini);
   double distance_neighbor=0.0;

   // Set as current best
   double best = distance;
        
   // Loop until system has cooled
   while (temp > 1) 
    {
     // Create new neighbour 3d model
     prot3d_model = prot3d_ini.clone();

     // Get a random conformation for this new neighbor
     prot3d_model.alterConformation();
            
     // Get energy of solutions
     currentEnergy = calcSimilarity(calcRCC(args[2]), calcRCC(prot3d_ini));;
     neighbourEnergy = calcSimilarity(calcRCC(args[2]), calcRCC(prot3d_model));;

     // Decide if we should accept the neighbour
     if (acceptanceProbability(currentEnergy, neighbourEnergy, temp) > Math.random()) {
       prot3d_ini = prot3d_model;
      }

     // Keep track of the best solution found
     if (neighbourEnergy < best) 
      {
       best = neighbourEnergy;
      }
            
     // Cool system
     temp *= 1-coolingRate;
    }

   System.out.println("Final solution distance: " + best.getDistance());
   System.out.println("Tour: " + best);
  }
}

