import java.util.HashMap;
import java.util.ArrayList;

abstract class AARotamers {

	protected static HashMap<String,Integer> sizeMap = new HashMap<String,Integer>();
	protected static HashMap<String,HashMap<Integer, ArrayList>> aminoMap  = new HashMap<String,HashMap<Integer, ArrayList>>();

	public static String getPDB(String aminoType, int rotNum, String chain, int aminoStart, int atomStart){
		int size = getSize(aminoType);
		HashMap<Integer, ArrayList> rotamers = getRotamers(aminoType);
		String pdb = "";
		String atNum,  amNum = "";
		String id, amino, x, y, z, type = "";
		int space = 0;
		for(int j = 0; j < size; j++){
			id = (String)((ArrayList)(rotamers.get(rotNum)).get(j)).get(0);
			amino = (String)((ArrayList)(rotamers.get(rotNum)).get(j)).get(1);
			x = (String)((ArrayList)(rotamers.get(rotNum)).get(j)).get(2);
			y = (String)((ArrayList)(rotamers.get(rotNum)).get(j)).get(3);
			z = (String)((ArrayList)(rotamers.get(rotNum)).get(j)).get(4);
			type = (String)((ArrayList)(rotamers.get(rotNum)).get(j)).get(5);
			atNum = String.valueOf(atomStart + j);
      space = 5 - atNum.length();
      for(int i=0; i < space; i++) atNum = " " + atNum;
			amNum = String.valueOf(aminoStart);
      space = 4 - amNum.length();
      for(int i=0; i < space; i++) amNum = " " + amNum;
					 //ATOM     10  CB  ALA A   2       1.438   0.894  -0.343  1.00 20.00           C
			pdb += "ATOM  " + atNum + id + amino + " " + chain + amNum + "    " + x + y + z + type;
			if(j < size -1){
 				pdb += "\n";
			}
		}
		return pdb;
	}

	public static int getSize(String aminoType){
		return sizeMap.get(aminoType);
	}

	protected static HashMap<Integer, ArrayList> getRotamers(String aminoType){
		return aminoMap.get(aminoType);
	}

	public static void main(String[] args) throws Exception {
		String am = args[0];
		System.out.println("main");
		System.out.println(aminoMap);
		System.out.println(getSize(am));
		System.out.println(getPDB(am, 1, "A", 2, 6));
	}
}
