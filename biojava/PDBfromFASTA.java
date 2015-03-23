package rccto3d;

import java.util.LinkedHashMap;
import java.io.File;
import java.io.IOException;

import org.biojava3.core.sequence.io.FastaReaderHelper;
import org.biojava3.core.sequence.ProteinSequence;


public class PDBfromFASTA {

	private String PDB;





	private void getPDB(ProteinSequence sequence){
		//for(){
		//}
	}
	
	private ProteinSequence readFasta(String filename, String identifier) throws Exception{
		LinkedHashMap<String, ProteinSequence> fasta = FastaReaderHelper.readFastaProteinSequence(new File(filename));
		return fasta.get(identifier);
	}

	private String getRotamer(String amino, int rotNum, String chain, int aminoStart, int atomStart){

		//return rotamers.ALARotamers.getPDB(rotNum, chain, aminoStart, atomStart);
		//switch(amino){
		//	case "A":
		//		return ALARotamers.getPDB(rotNum, chain, aminoStart, atomStart);
		//}
		return null;

	}

	public static void main(String[] args) throws Exception {
		LinkedHashMap<String, ProteinSequence> fasta = FastaReaderHelper.readFastaProteinSequence(new File(args[0]));
		System.out.println(fasta);
	}
}
