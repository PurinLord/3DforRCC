import java.util.LinkedHashMap;
import java.io.File;

import org.biojava3.core.sequence.io.FastaReaderHelper;
import org.biojava3.core.sequence.ProteinSequence;

public class PDBfromFASTA {






	public static void main(String[] args) throws Exception {
		LinkedHashMap<String, ProteinSequence> fasta = FastaReaderHelper.readFastaProteinSequence(new File(args[0]));
		System.out.println(fasta);
	}
}
