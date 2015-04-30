import java.util.HashMap;
import java.util.ArrayList;

public class AARotamers {

	private int size = 0;
	private HashMap<Integer,ArrayList> rotamers  = new HashMap<Integer,ArrayList>();
	private String aminoName = "";

	public String getPDB(int rotNum, String chain, int aminoStart, int atomStart){
		String pdb = "";
		String atNum,  amNum = "";
		String id, amino, x, y, z, type = "";
		int space = 0;
		for(int j = 0; j < size; j++){
			id = (String)((ArrayList)((ArrayList)rotamers.get(rotNum)).get(j)).get(0);
			amino = aminoName;
			x = (String)((ArrayList)((ArrayList)rotamers.get(rotNum)).get(j)).get(1);
			y = (String)((ArrayList)((ArrayList)rotamers.get(rotNum)).get(j)).get(2);
			z = (String)((ArrayList)((ArrayList)rotamers.get(rotNum)).get(j)).get(3);
			type = (String)((ArrayList)((ArrayList)rotamers.get(rotNum)).get(j)).get(4);
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

	public int getSize(){
		return this.size;
	}

	protected void setSize(int size){
		this.size = size;
	}

	protected void setRotamers(HashMap<Integer,ArrayList> rotamer){
		this.rotamers  = rotamer;
	}

	protected void setAminoName(String aminoName){
		this.aminoName = aminoName;
	}

	public String getAminoName(String aminoName){
		return this.aminoName;
	}

	public static void main(String[] args) throws Exception {
		String am = args[0];
		System.out.println("main");
		ALARotamers ala = new ALARotamers();
		System.out.println(ala.getSize());
		System.out.println(ala.getPDB(1, "A", 12, 63));
	}
}
