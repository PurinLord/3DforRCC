package rccto3d;

import java.util.LinkedHashMap;
import java.io.File;
import java.io.PrintWriter;
import java.io.IOException;

import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.Chain;
import org.biojava3.core.sequence.io.FastaReaderHelper;
import org.biojava3.core.sequence.ProteinSequence;

import org.biojava.bio.structure.AminoAcid;

import rccto3d.rotamers.*;


public class PDBfromFASTA {

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

	//TODO
	/*
	private int defaultASP = 1;
	private int defaultGLU = 1;
	private int defaultILE = 1;
	private int defaultPHE = 1;
	private int defaultTRP = 1;
	private int defaultALA = 1;
	private int defaultCYH = 1;
	private int defaultGLY = 1;
	private int defaultLEU = 1;
	private int defaultPRO = 1;
	private int defaultTYR = 1;
	private int defaultARG = 1;
	private int defaultCYS = 1;
	private int defaultHID = 1;
	private int defaultLYS = 1;
	private int defaultSER = 1;
	private int defaultVAL = 1;
	private int defaultASN = 1;
	private int defaultGLN = 1;
	private int defaultHIE = 1;
	private int defaultMET = 1;
	private int defaultTHR = 1;
	*/

	//public static final double ALFA_RIGHT_PHI = -0.9948; //-57°;
	//public static final double ALFA_RIGHT_PSI = -0.8203; //-47°;
	public static final double ALFA_RIGHT_PHI = -Math.PI*5.0/18.0; //-50°;
	public static final double ALFA_RIGHT_PSI = -Math.PI*5.0/18.0; //-50°;
	public static final double BETHA_ANTIPARALLEL_PHI = -2.426; //-139°;
	public static final double BETHA_ANTIPARALLEL_PSI = 2.356; //135°;
	public static final double NONE_PHI = 3.1426; //180°;
	public static final double NONE_PSI = 3.1426; //180°;

	public static void main(String[] args){
	try{
		PDBfromFASTA pff = new PDBfromFASTA();
		String pdb = pff.pdbFromFile(args[0], args[1]);
		//System.out.println(pdb);
		
		Trans.writePDB("datos/out.pdb", pff.shapeProtein(pdb, "alpha"));
		
		//Trans.writePDB("datos/tst9.pdb", pff.shapeProtein("tst.pdb", "alpha", true));
	}catch(Exception e){
  	e.printStackTrace();
	}
	}

	public Structure shapeProtein(String fileName, String shape, boolean keepFile) throws Exception{
		Structure struc = Trans.readPDB(fileName);
		Chain chain = struc.getChain(0);
		double phi = 0.0;
		double psi = 0.0;
			AminoAcid a1;
			AminoAcid a2;
		switch(shape){
			case "alpha":
				phi = ALFA_RIGHT_PHI;
				psi = ALFA_RIGHT_PSI;
				break;
			case "betha":
				phi = BETHA_ANTIPARALLEL_PHI;
				psi = BETHA_ANTIPARALLEL_PSI;
				break;
			case "none":
				phi = NONE_PHI;
				psi = NONE_PSI;
				break;
		}
		for(int i = 0; i < chain.getAtomLength() - 1; i++){
			Trans.makeBondTrans(chain, i);
			Trans.setPhi(chain, i, phi);
			Trans.setPsi(chain, i, psi);
				a1 = (AminoAcid)chain.getAtomGroup(i);
				a2 = (AminoAcid)chain.getAtomGroup(i + 1);
				System.out.println(Trans.getDihedral(a1, a2));
		}
		return struc;
	}

	public Structure shapeProtein(String pdb, String shape) throws Exception{
		PrintWriter writer = new PrintWriter("temp.pdb");
		writer.println(pdb);
		writer.close();
		return shapeProtein("temp.pdb", shape, true);
	}

	public String pdbFromFile(String fileName, String identifier, int from, int to, String chain, int startAmino, int startAtom) throws Exception{
		return pdbFromSequende(readFasta(fileName, identifier), from, to, chain, startAmino, startAtom);
	}

	public String pdbFromFile(String fileName, String identifier, String chain, int startAmino, int startAtom) throws Exception{
		return pdbFromSequende(readFasta(fileName, identifier), chain, startAmino, startAtom);
	}

	public String pdbFromFile(String fileName, String identifier) throws Exception{
		return pdbFromSequende(readFasta(fileName, identifier));
	}

	//starts at 1
	//ends at 0 or > length
	//start amino 1
	//start atom 1
	public String pdbFromSequende(ProteinSequence seq, int from, int to, String chain, int startAmino, int startAtom){
		String pdb = "";
		String amino;
		AARotamers curRotamer;
		int atomNum = startAtom;
		int aminoNum = startAmino - 1;
		int i;
		if (to > seq.getLength() || to == 0){
			to = seq.getLength();
		}
		for(i = 1; i <= to; i++){
				amino = seq.getCompoundAt(i + from - 1).toString();
				curRotamer = getRotamer(amino);
				pdb += curRotamer.getPDB(1, chain, i + aminoNum, atomNum) + "\n";
				atomNum += curRotamer.getSize();
		}
		return pdb;
	}

	public String pdbFromSequende(ProteinSequence seq, String chain, int startAmino, int startAtom){
		return pdbFromSequende(seq, 1, 0, chain, startAmino, startAtom);
	}

	public String pdbFromSequende(ProteinSequence seq){
		return pdbFromSequende(seq, "A", 1, 1);
	}

	//TODO
	public void setDefaultRotamer(String amino, int defRot){
		//String this.amino.getRotamers(amino);
	}
	
	private ProteinSequence readFasta(String filename, String identifier) throws Exception{
		LinkedHashMap<String, ProteinSequence> fasta = FastaReaderHelper.readFastaProteinSequence(new File(filename));
		return fasta.get(identifier);
	}

	//TODO
	private AARotamers getRotamer(String amino){
		switch(amino){
			case "A":
				return ALA;
			case "C":
				return CYS;
			case "D":
				return ASP;
			case "E":
				return GLU;
			case "F":
				return PHE;
			case "G":
				return GLY;
			case "H":
				return HIE;
			case "I":
				return ILE;
			case "K":
				return LYS;
			case "L":
				return LEU;
			case "M":
				return MET;
			case "N":
				return ASN;
			case "P":
				return PRO;
			case "Q":
				return GLN;
			case "R":
				return ARG;
			case "S":
				return SER;
			case "T":
				return THR;
			case "V":
				return VAL;
			case "W":
				return TRP;
			case "Y":
				return TYR;
		}
		System.out.println(">" + amino);
		return null;
	}
}
