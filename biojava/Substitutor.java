package rccto3d;

import java.io.IOException;
import java.lang.Integer;

import org.biojava.bio.structure.io.PDBFileReader;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.Chain;

import org.biojava.bio.structure.AminoAcid;

public class Substitutor{

Structure fitStruc = null;
PDBfromFASTA pff = null;

public Substitutor(Structure struc){
	fitStruc = struc;
}

public Structure fakeSubstitute(Structure outStruc, int startAt){
	//System.out.println(fitStruc.getChain(0).getAtomLength() + " " + outStruc.getChain(0).getAtomLength());
	for(int i = 0; i < outStruc.getChain(0).getAtomLength() -1; i ++){
		Chain chain = fitStruc.getChain(0);
		try{
			double angle = Trans.getPhiFix((AminoAcid)chain.getAtomGroup(i + startAt), (AminoAcid)chain.getAtomGroup(i + startAt + 1));
			Trans.setPhi(outStruc.getChain(0), i, angle);
		}catch(Exception e){
  		e.printStackTrace();
		}
	}
	return outStruc;
}

public Structure fakeSubstitute(String fileName, String identifier, int startAt){
	pff = new PDBfromFASTA();
	try{
		Structure outStruc = pff.makeProtein(pff.pdbFromFile(fileName, identifier));
	return fakeSubstitute(outStruc, startAt);
	}catch(Exception e){
  	e.printStackTrace();
		return null;
	}
}

public static void main(String args[]){
	Substitutor sub = new Substitutor(Trans.readPDB(args[0]));
	Structure out = sub.fakeSubstitute(Trans.readPDB(args[1]), Integer.parseInt(args[2]));
	Trans.writePDB("simDatos/out.pdb", out);
}

}
