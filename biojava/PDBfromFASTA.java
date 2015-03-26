package rccto3d;

import java.util.LinkedHashMap;
import java.io.File;
import java.io.IOException;

import org.biojava3.core.sequence.io.FastaReaderHelper;
import org.biojava3.core.sequence.ProteinSequence;

import rccto3d.rotamers.*;


public class PDBfromFASTA {

	private String pdb = "";

	private final ASPRotamers ASP = new ASPRotamers();
	private final GLURotamers GLU = new GLURotamers();
	private final ILERotamers ILE = new ILERotamers();
	private final PHERotamers PHE = new PHERotamers();
	private final TRPRotamers TRP = new TRPRotamers();
	private final ALARotamers ALA = new ALARotamers();
	private final CYHRotamers CYH = new CYHRotamers();
	private final GLYRotamers GLY = new GLYRotamers();
	private final LEURotamers LEU = new LEURotamers();
	private final PRORotamers PRO = new PRORotamers();
	private final TYRRotamers TYR = new TYRRotamers();
	private final ARGRotamers ARG = new ARGRotamers();
	private final CYSRotamers CYS = new CYSRotamers();
	private final HIDRotamers HID = new HIDRotamers();
	private final LYSRotamers LYS = new LYSRotamers();
	private final SERRotamers SER = new SERRotamers();
	private final VALRotamers VAL = new VALRotamers();
	private final ASNRotamers ASN = new ASNRotamers();
	private final GLNRotamers GLN = new GLNRotamers();
	private final HIERotamers HIE = new HIERotamers();
	private final METRotamers MET = new METRotamers();
	private final THRRotamers THR = new THRRotamers();

	public String makeSequence(String fileName, String identifier) throws Exception{
		getPDB(readFasta(fileName, identifier));
		return pdb;
	}

	private void getPDB(ProteinSequence seq){
		String amino;
		int atomNum = 0;
		int i;
		for(i = 1; i <= seq.getLength(); i++){
				amino = seq.getCompoundAt(i).toString();
				pdb += getRotamer(amino, 1, "A", i, atomNum) + "\n";
				//atomNum += 
		}
	}
	
	private ProteinSequence readFasta(String filename, String identifier) throws Exception{
		LinkedHashMap<String, ProteinSequence> fasta = FastaReaderHelper.readFastaProteinSequence(new File(filename));
		return fasta.get(identifier);
	}

	private String getRotamer(String amino, int rotNum, String chain, int aminoStart, int atomStart){

		switch(amino){
			case "A":
				return ALA.getPDB(rotNum, chain, aminoStart, atomStart);
			case "C":
				return CYS.getPDB(rotNum, chain, aminoStart, atomStart);
			case "D":
				return ASP.getPDB(rotNum, chain, aminoStart, atomStart);
			case "E":
				return GLU.getPDB(rotNum, chain, aminoStart, atomStart);
			case "F":
				return PHE.getPDB(rotNum, chain, aminoStart, atomStart);
			case "G":
				return GLY.getPDB(rotNum, chain, aminoStart, atomStart);
			case "H":
				return HIE.getPDB(rotNum, chain, aminoStart, atomStart);
			case "I":
				return ILE.getPDB(rotNum, chain, aminoStart, atomStart);
			case "K":
				return LYS.getPDB(rotNum, chain, aminoStart, atomStart);
			case "L":
				return LEU.getPDB(rotNum, chain, aminoStart, atomStart);
			case "M":
				return MET.getPDB(rotNum, chain, aminoStart, atomStart);
			case "N":
				return ASN.getPDB(rotNum, chain, aminoStart, atomStart);
			case "P":
				return PRO.getPDB(rotNum, chain, aminoStart, atomStart);
			case "Q":
				return GLN.getPDB(rotNum, chain, aminoStart, atomStart);
			case "R":
				return ARG.getPDB(rotNum, chain, aminoStart, atomStart);
			case "S":
				return SER.getPDB(rotNum, chain, aminoStart, atomStart);
			case "T":
				return THR.getPDB(rotNum, chain, aminoStart, atomStart);
			case "V":
				return VAL.getPDB(rotNum, chain, aminoStart, atomStart);
			case "W":
				return TRP.getPDB(rotNum, chain, aminoStart, atomStart);
			case "Y":
				return TYR.getPDB(rotNum, chain, aminoStart, atomStart);
		}
		System.out.println("[]" + amino);
		return null;
	}

	public static void main(String[] args){
	try{
		PDBfromFASTA pff = new PDBfromFASTA();
		System.out.println(pff.makeSequence(args[0], args[1]));
	}catch(Exception e){
  	e.printStackTrace();
	}
	}
}
