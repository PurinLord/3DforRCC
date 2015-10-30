package rccto3d;

import java.util.LinkedHashMap;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Random;

import org.biojava.bio.structure.io.PDBFileReader;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.Chain;
import org.biojava3.core.sequence.io.FastaReaderHelper;
import org.biojava3.core.sequence.ProteinSequence;

import org.biojava.bio.structure.AminoAcid;

import rccto3d.rotamers.*;


public class PDBfromFASTA {

	private int defaultASP = 1;
	private int defaultGLU = 1;
	private int defaultILE = 1;
	private int defaultPHE = 1;
	private int defaultTRP = 1;
	private int defaultALA = 1;
	//private int defaultCYH 1;;
	private int defaultGLY = 1;
	private int defaultLEU = 1;
	private int defaultPRO = 1;
	private int defaultTYR = 1;
	private int defaultARG = 1;
	private int defaultCYS = 1;
	//private int defaultHID 1;
	private int defaultLYS = 1;
	private int defaultSER = 1;
	private int defaultVAL = 1;
	private int defaultASN = 1;
	private int defaultGLN = 1;
	private int defaultHIS = 1;
	private int defaultMET = 1;
	private int defaultTHR = 1;

	private final ASPRotamers ASP = new ASPRotamers(defaultASP);
	private final GLURotamers GLU = new GLURotamers(defaultGLU);
	private final ILERotamers ILE = new ILERotamers(defaultILE);
	private final PHERotamers PHE = new PHERotamers(defaultPHE);
	private final TRPRotamers TRP = new TRPRotamers(defaultTRP);
	private final ALARotamers ALA = new ALARotamers(defaultALA);
	//private final CYHRotamers CYH = new CYHRotamers();
	private final GLYRotamers GLY = new GLYRotamers(defaultGLY);
	private final LEURotamers LEU = new LEURotamers(defaultLEU);
	private final PRORotamers PRO = new PRORotamers(defaultPRO);
	private final TYRRotamers TYR = new TYRRotamers(defaultTYR);
	private final ARGRotamers ARG = new ARGRotamers(defaultARG);
	private final CYSRotamers CYS = new CYSRotamers(defaultCYS);
	//private final HIDRotamers HID = new HIDRotamers();
	private final LYSRotamers LYS = new LYSRotamers(defaultLYS);
	private final SERRotamers SER = new SERRotamers(defaultSER);
	private final VALRotamers VAL = new VALRotamers(defaultVAL);
	private final ASNRotamers ASN = new ASNRotamers(defaultASN);
	private final GLNRotamers GLN = new GLNRotamers(defaultGLN);
	private final HISRotamers HIS = new HISRotamers(defaultHIS);
	private final METRotamers MET = new METRotamers(defaultMET);
	private final THRRotamers THR = new THRRotamers(defaultTHR);

	//public static final double ALFA_RIGHT_PHI = -0.9948; //-57°;
	//public static final double ALFA_RIGHT_PSI = -0.8203; //-47°;
	public static final double ALFA_RIGHT_PHI = -50; //-Math.PI*5.0/18.0;
	public static final double ALFA_RIGHT_PSI = -50; //-Math.PI*5.0/18.0;
	public static final double BETHA_ANTIPARALLEL_PHI = -139; //-2.426;
	public static final double BETHA_ANTIPARALLEL_PSI =  135; // 2.356;
	public static final double NONE_PHI = 180; //3.1426;
	public static final double NONE_PSI = 180; //3.1426;
	private String tempDirectory;

	public static void main(String[] args){
		try{
			PDBfromFASTA pff = new PDBfromFASTA();
			String pdb = pff.pdbFromFile(args[0], args[1]);
			//System.out.println(pdb);
			
			writePDB("datos/out.pdb", pff.shapeProtein(pdb, args[2], false));
			
			//Trans.writePDB("datos/tst9.pdb", pff.shapeProtein("tst.pdb", "alpha", true));
		}catch(Exception e){
  		e.printStackTrace();
		}
	}

	//Todavía no estoy seguro de porqué hice esta clase TAN confusa, todo llama a todo por todos lados O.o

	public PDBfromFASTA (){
		tempDirectory = "temp.pdb"; 
		GLY.setRotNum(defaultGLY);
	}

	public PDBfromFASTA (String directory){
		tempDirectory = directory + "temp.pdb"; 
	}

	/** First it puts together the aminoacids usin makeProtein, then it folds it in to the chosen shape
	*/
	public Structure shapeProtein(String fileName, String shape, boolean keepFile) throws Exception{
		Structure struc = makeProtein(fileName, keepFile);
		Chain chain = struc.getChain(0);
		double phi = -20.0;
		double psi = -20.0;
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
			//System.out.println("phi " + phi + " psi " + psi);
		for(int i = 0; i < chain.getAtomLength() - 1; i++){
			Trans.setPhi(chain, i, phi);
			Trans.setPsi(chain, i, psi);
		}
		return struc;
	}

	public Structure shapeProtein(String pdb) throws Exception{
		File file = new File(tempDirectory);
		if (!file.exists()) {
			file.createNewFile();
		}
		synchronized(file) {
			BufferedWriter bw = new BufferedWriter(new FileWriter(file.getAbsoluteFile()));
			bw.write(pdb);
			bw.close();
		}
		return shapeProtein(tempDirectory, "none", true);
	}

	public Structure randomShapeProtein(String pdb, long seed) throws Exception{
		Structure struc = makeProtein(pdb);
		Chain chain = struc.getChain(0);
		Random rdm = new Random(seed);
		for(int i = 0; i < chain.getAtomLength() - 1; i++){
			Trans.setPhi(chain, i, rdm.nextDouble() * 360);
			Trans.setPsi(chain, i, rdm.nextDouble() * 360);
		}
		return struc;
	}

	/** Pusts together all the amioacids in the pdb with the correct bond distance and omega angle
	 */
	public Structure makeProtein(String fileName, boolean keepFile) throws Exception{
		Structure struc = readPDB(fileName);
		Chain chain = struc.getChain(0);
		for(int i = 0; i < chain.getAtomLength() - 1; i++){
			Trans.makeBondTrans(chain, i);
		}
		if(!keepFile){
			try{
				File file = new File(fileName);
    		file.delete();
			}catch(Exception e){
    		e.printStackTrace();
    	}
		}
		return struc;
	}

	public Structure makeProtein(String pdb) throws Exception{
		File file = new File(tempDirectory);
		if (!file.exists()) {
			file.createNewFile();
		}
		synchronized(file) {
			BufferedWriter bw = new BufferedWriter(new FileWriter(file.getAbsoluteFile()));
			bw.write(pdb);
			bw.close();
		}
		return makeProtein(tempDirectory, false);
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
				return HIS;
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
		//System.out.println(">" + amino);
		return null;
	}

	public static void writePDB(String filename, Structure structure){
		try{
			File file = new File(filename);
			if (!file.exists()) {
				file.createNewFile();
			}
			synchronized(file) {
				BufferedWriter bw = new BufferedWriter(new FileWriter(file.getAbsoluteFile()));
				bw.write(structure.toPDB().trim());
				bw.close();
			}
		}catch(Exception e){
			e.printStackTrace();
		}
	}

	public static Structure readPDB(String filename){
 		PDBFileReader pdbreader = new PDBFileReader();
		Structure structure = null;
		File file = new File(filename);
    try{
			synchronized(file) {
    		structure = pdbreader.getStructure(file);
			}
    } catch (Exception e) {
    	e.printStackTrace();
    }
    return structure;
	}

	//TODO
	public void setDefaultRotamer(String amino, int defRot){
		//String this.amino.getRotamers(amino);
	}
}
