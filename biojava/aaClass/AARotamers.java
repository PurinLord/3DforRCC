import java.util.HashMap;
import java.util.ArrayList;

abstract class AARotamers {

	protected static int size = 0;
	protected static HashMap<Integer,ArrayList> rotamers  = new HashMap<Integer,ArrayList>();

	public static String getPDB(int rotNum, String chain, int aminoStart, int atomStart){
		String pdb = "";
		String atNum,  amNum = "";
		String id, amino, x, y, z, type = "";
		int space = 0;
		for(int j = 0; j < size; j++){
			id = (String)((ArrayList)((ArrayList)rotamers.get(rotNum)).get(j)).get(0);
			amino = (String)((ArrayList)((ArrayList)rotamers.get(rotNum)).get(j)).get(1);
			x = (String)((ArrayList)((ArrayList)rotamers.get(rotNum)).get(j)).get(2);
			y = (String)((ArrayList)((ArrayList)rotamers.get(rotNum)).get(j)).get(3);
			z = (String)((ArrayList)((ArrayList)rotamers.get(rotNum)).get(j)).get(4);
			type = (String)((ArrayList)((ArrayList)rotamers.get(rotNum)).get(j)).get(5);
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

	public static int getSize(){
		return size;
	}

	protected static void setSize(int newSize){
		size = newSize;
	}

	protected static void setRotamers(HashMap<Integer,ArrayList> newRotamer){
		rotamers  = newRotamer;
	}

	public static void main(String[] args) throws Exception {
		System.out.println("main");
		System.out.println(rotamers);
		System.out.println(getSize());
		System.out.println(getPDB(1, "A", 2, 6));
	}
}
