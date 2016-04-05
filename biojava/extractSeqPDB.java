import java.io.*;
import java.util.*;

import java.io.File;
import java.io.BufferedWriter;
import java.lang.StringBuilder;

public class extractSeqPDB
{
	String targetDir;
	public extractSeqPDB(String targetDir){
		this.targetDir = targetDir;
	}

 public static String translate(String aa)
 {
  String result="";
  if(aa.equalsIgnoreCase("ALA")) result="A";
  if(aa.equalsIgnoreCase("CYS")) result="C";
  if(aa.equalsIgnoreCase("ASP")) result="D";
  if(aa.equalsIgnoreCase("GLU")) result="E";
  if(aa.equalsIgnoreCase("PHE")) result="F";
  if(aa.equalsIgnoreCase("GLY")) result="G";
  if(aa.equalsIgnoreCase("HIS")) result="H";
  if(aa.equalsIgnoreCase("ILE")) result="I";
  if(aa.equalsIgnoreCase("LYS")) result="K";
  if(aa.equalsIgnoreCase("LEU")) result="L";
  if(aa.equalsIgnoreCase("MET")) result="M";
  if(aa.equalsIgnoreCase("ASN")) result="N";
  if(aa.equalsIgnoreCase("PRO")) result="P";
  if(aa.equalsIgnoreCase("GLN")) result="Q";
  if(aa.equalsIgnoreCase("ARG")) result="R";
  if(aa.equalsIgnoreCase("SER")) result="S";
  if(aa.equalsIgnoreCase("THR")) result="T";
  if(aa.equalsIgnoreCase("VAL")) result="V";
  if(aa.equalsIgnoreCase("TRP")) result="W";
  if(aa.equalsIgnoreCase("TYR")) result="Y";
  return result;
 } // end of method translate

 public void extract(String inPDB, String chain) throws IOException, FileNotFoundException{
    //String file = inPDB;
		//String chain = args[1]
		String line="", residue="", currentChain="", number="";
    StringBuffer print = new StringBuffer();
    StringTokenizer st = null;
    BufferedReader infile = new BufferedReader(new FileReader(inPDB));
    Vector seq = new Vector();
		String fasta = ">a\n";
		String pdb = "";
		String name = inPDB.substring(inPDB.lastIndexOf("/"), inPDB.length()-4);

    while((line=infile.readLine())!=null){
      if(line.startsWith("ATOM")){
        st = new StringTokenizer(line);
        st.nextElement(); st.nextElement(); 
		  String checkRes = (String)st.nextElement();
		  if(checkRes.length() > 3){
			  residue = checkRes.substring(3);
		  }else{
        		residue=(String)st.nextElement();
		  }

        currentChain=(String)st.nextElement();

        if(Character.isDigit(currentChain.charAt(0)) && chain.equalsIgnoreCase("none")){
          if(!seq.contains(currentChain)){
            seq.add(currentChain);
            print.append(translate(residue));
            if(print.length()==50){
              //System.out.println(print);
							fasta += print + "\n";
              print.delete(0,50);
             }
           }
         }
        if(Character.isLetter(currentChain.charAt(0)) && currentChain.equalsIgnoreCase(chain)){
          number=(String)st.nextElement();
			 String subAmino = "";
			String resID = "";
			if(residue.length() == 3){
				resID = translate(residue);
				line += "\n";
			}else{
				subAmino = residue.substring(0,1);
				if(subAmino.equals("A")){
					resID = translate(residue.substring(1));
					StringBuilder edit = new StringBuilder(line);
					edit.setCharAt(16, ' ');
					line = edit.toString() + "\n";
				}else{
					line = "";
				}
			 }
				pdb += line;

          if(!seq.contains(number)){
            seq.add(number);
            print.append(resID);
            if(print.length()==50){
              //System.out.println(print.toString());
					fasta += print + "\n";
              print.delete(0,50); 
				}
			 }
		  }
        while(st.hasMoreElements()) st.nextElement();
		}
	 }
    infile.close();
    //System.out.println(print.toString());
		try{
			File f = new File(targetDir);
			if (!f.exists()){
				f.mkdirs();
			}
			File file = new File(targetDir + name +".fa");
			BufferedWriter bw = new BufferedWriter(new FileWriter(file.getAbsoluteFile()));
			fasta += print;
			bw.write(fasta.trim());
			bw.close();

			file = new File(targetDir + name +".pdb");
			bw = new BufferedWriter(new FileWriter(file.getAbsoluteFile()));
			bw.write(pdb.trim());
			bw.close();
		}catch(Exception e){
			e.printStackTrace();
		}
    //System.out.println(fasta);
    //System.out.println(pdb);
 }

 public static void main(String[] args)
 {
  if(args.length==3)
   {
		 try{
			extractSeqPDB extr = new extractSeqPDB(args[2]);
			extr.extract(args[0], args[1]);
		 }catch(Exception e){
    	e.printStackTrace();
		 }
   }
  else
   {
    System.err.println("Usage:\njava extractSeqPDB <pdb> <chain> <targetDir>");
    System.err.println("<pdb> pdb file name");
    System.err.println("<chain> an String referring to a chain in <pdb>.");
    System.err.println("The program will extract the sequence from the <chain> specified");
    System.err.println("in <pdb> the sequence. For that, the program will read the");
    System.err.println("artomic coordinates and translate the 3-letter amino acid code");
    System.err.println("used in PDB to 1-letter code.");
    System.err.println("The ourput will be sent to the standard output.");
    System.err.println("");
   }
 }
}
