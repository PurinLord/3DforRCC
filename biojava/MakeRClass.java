import java.io.*;
import java.util.Arrays;

public class MakeRClass{

public static void main(String[] args) throws Exception {
  if(args.length==1){

	PrintWriter writer;
	String amino, rotNum;
	String fileName;
  String line, ID, atom, x, y, z, type;
	int size;
  line = ID = atom = x = y = z = type = fileName = amino = rotNum = "";
	String prog = "";
	size = 0;

	String head = "import java.util.HashMap;\n" +
								"import java.util.ArrayList;\n\n" +
								"public class ";
	String headEnd =  " extends AARotamers {\n" +
										"\tprivate static HashMap<Integer,ArrayList> rotamers  = new HashMap<Integer,ArrayList>();\n" +
										"\tprivate static ArrayList<ArrayList> aminoAcid;\n" +
										"\tprivate static ArrayList<String> atom;\n\n" +
										"\tstatic {\n" +
										"\t\tsizeMap.put(\"YYY\", XX);\n" +
										"\t\taminoMap.put(\"YYY\", rotamers);\n";

	String startAminoAcid = "\t\taminoAcid = new ArrayList<ArrayList>();\n";
	String startAtom = 		"\t\tatom = new ArrayList<String>();\n";
	String putAminoAcid = 	"\t\trotamers.put(";
	String putAminoAcidEnd =",aminoAcid);\n";
	String addAtom = 			"\t\taminoAcid.add(atom);\n";
	String addPar = 			"atom.add(";
	String addParEnd = 		"); ";
	String nline = 				"\n";

	File directory = new File(args[0]);
	File[] fList = directory.listFiles();
	Arrays.sort(fList);

	for (File file : fList){
  	ID = atom = x = y = z = type = "";
		BufferedReader infile = new BufferedReader(new FileReader(file));
		fileName = file.getName();
		if (!amino.equals(fileName.substring(0,3))){
			if(!amino.equals("")){
				prog += "\t}\n";
				prog += "}\n";
<<<<<<< HEAD
				prog = prog.replaceFirst("XX", String.valueOf(size));
				prog = prog.replaceAll("YYY", amino);
				writer = new PrintWriter(dirName + amino + "Rotamers.java", "UTF-8");
=======
				prog = prog.replaceFirst("XXX", String.valueOf(size));
				writer = new PrintWriter("aaClass/" + amino + "Rotamers.java", "UTF-8");
>>>>>>> parent of a4698c0... added package and minor changes
				writer.println(prog);
				writer.close();
				prog = "";
			}
			amino = fileName.substring(0,3);
			prog += head + amino + "Rotamers" + headEnd;
		}
		if(fileName.substring(10,11).equals(".")){
			rotNum = fileName.substring(9,10);
		}else{
		rotNum = fileName.substring(9,11);
		}
		prog += nline;
		prog += "\t\t//Rostamer " + rotNum + "\n";
		prog += startAminoAcid;
		prog += putAminoAcid + rotNum + putAminoAcidEnd;

		size = 0;
		while((line=infile.readLine())!=null){
			size ++;
			ID = line.substring(11,17);
			atom = line.substring(17,20);
			x = line.substring(30,38);
			y = line.substring(38,46);
			z = line.substring(46,54);
			type = line.substring(54,78);
			prog += nline;
			prog += startAtom;
			prog += "\t\t" + addPar + "\"" + ID + "\"" + addParEnd +
							addPar + "\"" + atom + "\"" + addParEnd +
							addPar + "\"" + x + "\"" + addParEnd +
							addPar + "\"" + y + "\"" + addParEnd +
							addPar + "\"" + z + "\"" + addParEnd +
							addPar + "\"" + type + "\"" + addParEnd + "\n"; 
			prog += addAtom;
		}
	}
	prog += "\t}\n";
	prog += "}\n";
<<<<<<< HEAD
	prog = prog.replaceFirst("XX", String.valueOf(size));
	prog = prog.replaceAll("YYY", amino);
	writer = new PrintWriter(dirName + amino + "Rotamers.java", "UTF-8");
=======
	prog = prog.replaceFirst("XXX", String.valueOf(size));
	writer = new PrintWriter("aaClass/" + amino + "Rotamers.java", "UTF-8");
>>>>>>> parent of a4698c0... added package and minor changes
	writer.println(prog);
	writer.close();
	//System.out.println(ID + "|" + atom + "|");
	//System.out.println(x + "|" + y + "|" + z + "|");
	}else {
		System.err.println("Usage:\njava MakeRClass <directory>");
	}
}
}
