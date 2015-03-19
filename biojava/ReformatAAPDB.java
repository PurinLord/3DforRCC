import java.io.*;
import java.util.*;

public class ReformatAAPDB
{
 public static void main(String[] args) throws Exception
  {
   if(args.length==1)
    {

		 File directory = new File(args[0]);
		 File[] fList = directory.listFiles();
		 for (File file : fList){
		 BufferedReader infile = new BufferedReader(new FileReader(file));
     //BufferedReader infile = new BufferedReader(new FileReader(args[0]));
     int i=0, emptySpaces=0, serialNo=1, resNum=1;
     String atom="ATOM  ";
     String line="", newLine="", fullLine="", serial="", atomName="", resName="", chainID="A", resSeq="   1", token="";

     while((line=infile.readLine())!=null)
      {
       if(line.startsWith("ATOM"))
        {
         atomName=line.substring(11,15);
         if(atomName.indexOf("H")==-1)
          {
           token=line.substring(27);
           newLine=atom;

           serial=String.valueOf(serialNo); serialNo=serialNo+1;
           emptySpaces=5-serial.length();
           for(i=0; i<emptySpaces; i++) serial=" "+serial;
           newLine=newLine+serial+" ";
          
           emptySpaces=4-atomName.length();
           for(i=0; i<emptySpaces; i++) atomName=" "+atomName;
					 //System.out.println(atomName + "|");
           newLine=newLine+atomName;

           resName=line.substring(17,20);
           emptySpaces=4-resName.length();
           for(i=0; i<emptySpaces; i++) resName=" "+resName;
           newLine=newLine+resName+" ";

           newLine=newLine+chainID;

           newLine=newLine+resSeq+" ";

           newLine=newLine+token;

           System.out.println(newLine);
					 fullLine += newLine + "\n";
          }
        }
      }
		 fullLine = fullLine.trim();
     System.out.println("ATOM     15  OE1 GLN A   2       5.366   6.856  36.946  1.00 31.80           O");
     //System.out.println(file.getPath());
		 PrintWriter writer = new PrintWriter(file.getPath(), "UTF-8");
		 writer.println(fullLine);
		 writer.close();
     infile.close();
    }
		}
   else
    {
     System.err.println("Usage:\njava ReformatAAPDB <infile>");
     System.err.println("<infile> file with PDB coordinates obtained from DYNAMEOMICS db; e.g., ARG_2014_1.pdb.txt.");
     System.err.println("The program will print out to the standard output the section corresponding");
     System.err.println("to the ATOm coordinates, whithout Hydrogens, with each residue being part of chain A.");
     System.err.println("");
     System.err.println("");
    }
  }
}
