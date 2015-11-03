	import java.io.File;
	import java.io.IOException;
	import java.io.PrintWriter;
	import java.net.URL;
	import java.lang.Math;
	
	//import org.biojava.bio.structure.*;
	import org.biojava.bio.structure.Structure;
	import org.biojava.bio.structure.Chain;
	import org.biojava.bio.structure.Group;
	import org.biojava.bio.structure.AminoAcid;
	import org.biojava.bio.structure.AminoAcidImpl;
	import org.biojava.bio.structure.Atom; import org.biojava.bio.structure.AtomImpl;
	import org.biojava.bio.structure.io.PDBFileReader;
	import org.biojava.bio.structure.Calc;
	import org.biojava.bio.structure.jama.Matrix;
	
	import java.io.*;
	import java.util.HashMap;
	import java.util.LinkedHashMap;
	import java.util.ArrayList;
	
	import org.biojava3.core.sequence.io.FastaReader;
	import org.biojava3.core.sequence.io.FastaReaderHelper;
	import org.biojava3.core.sequence.io.*;
	import org.biojava3.core.sequence.ProteinSequence;
	import org.biojava3.core.sequence.compound.AminoAcidCompound;
	import org.biojava3.core.sequence.compound.AminoAcidCompoundSet;
	import org.biojava3.core.sequence.compound.DNACompoundSet;
	import org.biojava3.core.sequence.compound.NucleotideCompound;
	import org.biojava3.core.sequence.DNASequence;
	//import org.biojava3.core.sequence.io.template
	//import org.biojava3.core.sequence.io.template.SequenceCreatorInterface;
	//import org.biojava3.core.sequence.io.template.SequenceHeaderParserInterface;
	
	//import org.biojava.bio.structure.Mutator;
	//import org.biojava.bio.seq.*;
	//import org.biojava.bio.symbol.*;
	
	import rccto3d.optimisation.SimulatedAnnealing3DProtFromRCC;
	import rccto3d.PDBfromFASTA;
	import rccto3d.optimisation.*;
	import rccto3d.rotamers.*;
	import rccto3d.*;

/** Makes a solid body rotatios of a spesific residue to change its dihedral angle
 */
public class TestRot {
  	static Structure struc;
  	static Structure struc2;
		static Chain chain;
		static Chain chain2;
		static double phi;
		static double psi;
		static double pi = 3.1416;
		static double trans = pi/180.0;
		static int res = 45;

		static double[] axis;
		static double[] axis2;
		static double angle;
		static double angle2;
		static double angle3;

		static AminoAcid a1;
		static AminoAcid a2;

		static final double A_PSI = -50;
		static final double A_PHI = -60;

  public static void main(String[] args) {

    try {

			angle = 180.0;
			//angle2 = -Math.PI*5.0/18.0;
			//angle2 = -50;
			angle2 = 0;

			struc = rccto3d.Trans.readPDB(args[0]);
			//struc = rccto3d.Trans.readPDB("datos/tst0.pdb");
			chain = struc.getChain(0);

			rccto3d.PDBfromFASTA pff = new rccto3d.PDBfromFASTA();

			SimulatedAnnealing3DProtFromRCC simA = new SimulatedAnnealing3DProtFromRCC("", "");
			//simA.calcInitialTemp(3, struc);


			//chain2 = (Chain)chain.clone();
			//		a1 = (AminoAcid)chain.getAtomGroup(0);
			//		System.out.println(a1.toString());
			//		System.out.println(a1.getAtom(0).toString());
			//		System.out.println(a1.getAtom(0).getPDBserial());
			//		System.out.println(a1.getResidueNumber().getSeqNum());

			//struc2 = rccto3d.Trans.readPDB(args[0]);

			//for(int i = 0; i < chain.getAtomLength() - 1; i++){
			//	a1 = (AminoAcid)chain.getAtomGroup(i);
			//	a2 = (AminoAcid)chain.getAtomGroup(i + 1);
			//		//System.out.println(a1.toString());
			//	angle = rccto3d.Trans.getPhi(a1, a2);
			//	angle2 = rccto3d.Trans.getPsi(a1, a2);
			//		System.out.println(angle + "  " + angle2);
			//		//System.out.println(rccto3d.Trans.getDihedral(a1, a2));
			//	//rccto3d.Trans.setOmega(chain, i, angle);
			//	//angle = Math.round(angle);
			//	//angle2 = Math.round(angle2);
			//		//System.out.println(angle + "  " + angle2);
			//	rccto3d.Trans.setPhi(chain2, i, angle);
			//	rccto3d.Trans.setPsi(chain2, i, angle2);
			//	a1 = (AminoAcid)chain2.getAtomGroup(i);
			//	a2 = (AminoAcid)chain2.getAtomGroup(i + 1);
			//	angle = rccto3d.Trans.getPhi(a1, a2);
			//	angle2 = rccto3d.Trans.getPsi(a1, a2);
			//	System.out.println(angle + "  " + angle2);
			//		//System.out.println(rccto3d.Trans.getDihedral(a1, a2));
			//	//rccto3d.Trans.makeBondTrans(chain, i);
			//		//System.out.println(rccto3d.Trans.getDihedral(a1, a2));
			//}
			//rccto3d.Trans.writePDB("datos/tst1.pdb", chain2.getParent());

			//////	Redondea ángulos
			//struc2 = (Structure)struc.clone();
			//for(int i = 0; i < chain.getAtomLength() - 1; i++){
			//	a1 = (AminoAcid)chain.getAtomGroup(i);
			//	a2 = (AminoAcid)chain.getAtomGroup(i + 1);
			//	angle = rccto3d.Trans.getPhi(a1, a2);
			//	angle2 = rccto3d.Trans.getPsi(a1, a2);
			//	System.out.println(angle + "  " + angle2);
			//	angle = Math.round(angle);
			//	angle2 = Math.round(angle2);
			//	rccto3d.Trans.setPhi(chain, i, angle);
			//	rccto3d.Trans.setPsi(chain, i, angle2);
			//	angle = rccto3d.Trans.getPhi(a1, a2);
			//	angle2 = rccto3d.Trans.getPsi(a1, a2);
			//	System.out.println(angle + "  " + angle2);
			//}
			//double dist[] = simA.calcEnergy(struc, struc2);
			//System.out.println(dist[0] + "  " + dist[1]);
			//rccto3d.PDBfromFASTA.writePDB("datos/tst1.pdb", chain.getParent());

			//////	Checa diferencia entre estructura original y la generada
			String pdb = pff.pdbFromFile(args[1], "a");
			Structure struc2 = pff.shapeProtein(pdb);
			Chain chain2 = struc2.getChain(0);
			//rccto3d.Trans.writePDB("simDatos/deltaEst_pre.pdb", chain2.getParent());
			for(int i = 0; i < chain.getAtomLength() - 1; i++){
			rccto3d.Trans.writePDB("datos/d_"+i+".pdb", chain2.getParent());
				a1 = (AminoAcid)chain.getAtomGroup(i);
				a2 = (AminoAcid)chain.getAtomGroup(i + 1);
				angle = rccto3d.Trans.getPhi(a1, a2);
				angle2 = rccto3d.Trans.getPsi(a1, a2);
				angle3 = rccto3d.Trans.getOmega(a1, a2);
				System.out.println(angle + "  " + angle2 + "  " + angle3);
				rccto3d.Trans.setPhi(chain2, i, angle);
				rccto3d.Trans.setPsi(chain2, i, angle2);
				//rccto3d.Trans.setOmega(chain2, i, angle3);
				a1 = (AminoAcid)chain2.getAtomGroup(i);
				a2 = (AminoAcid)chain2.getAtomGroup(i + 1);
				angle = rccto3d.Trans.getPhi(a1, a2);
				angle2 = rccto3d.Trans.getPsi(a1, a2);
				angle3 = rccto3d.Trans.getOmega(a1, a2);
				System.out.println(angle + "  " + angle2 + "  " + angle3);
			}
			rccto3d.Trans.writePDB("simDatos/deltaEst.pdb", chain2.getParent());
			System.out.println("> "+simA.calcEnergy(chain.getParent(),chain2.getParent())[0]);
			System.out.println("> "+simA.calcEnergy(chain.getParent(),chain2.getParent())[1]);

			////////	Checa diferencia entre omega
			//String pdb = pff.pdbFromFile(args[1], "a");
			//Structure struc2 = pff.shapeProtein(pdb);
			//Chain chain2 = struc2.getChain(0);
			//for(int i = 0; i < chain.getAtomLength() - 1; i++){
			//	a1 = (AminoAcid)chain.getAtomGroup(i);
			//	a2 = (AminoAcid)chain.getAtomGroup(i + 1);
			//	angle = rccto3d.Trans.getPhi(a1, a2);
			//	angle2 = rccto3d.Trans.getPsi(a1, a2);
			//	angle3 = rccto3d.Trans.getOmega(a1, a2);
			//	rccto3d.Trans.setPhi(chain2, i, angle);
			//	rccto3d.Trans.setPsi(chain2, i, angle2);
			//	a1 = (AminoAcid)chain2.getAtomGroup(i);
			//	a2 = (AminoAcid)chain2.getAtomGroup(i + 1);
			//	angle2 =rccto3d.Trans.getOmega(a1, a2);
			//	if(angle3>angle2){
			//	angle3 = angle3 - angle2; 
			//	}else{
			//	angle3 =angle2 - angle3;
			//	}
			//	if(angle3>180){
			//		angle3 = 360 - angle3;
			//	}
			//	System.out.println(angle3);
			//}

			////////	Regresa estructura a tripa
			//for(int i = 0; i < chain.getAtomLength() - 1; i++){
			//	rccto3d.Trans.setPhi(chain, i, 180);
			//	rccto3d.Trans.setPsi(chain, i, 180);
			//	a1 = (AminoAcid)chain.getAtomGroup(i);
			//	a2 = (AminoAcid)chain.getAtomGroup(i + 1);
			//	angle = rccto3d.Trans.getPhi(a1, a2);
			//	angle2 = rccto3d.Trans.getPsi(a1, a2);
			//	System.out.println(angle + "  " + angle2);
			//}
			//rccto3d.PDBfromFASTA.writePDB("simDatos/test_oxigen.pdb", chain.getParent());

			//for(int i = 0; i < chain.getAtomLength() - 1; i++){
			//	rccto3d.Trans.setPhi(chain, i, angle2);
			//	rccto3d.Trans.setPsi(chain, i, angle2);
			//		a1 = (AminoAcid)chain.getAtomGroup(i);
			//		a2 = (AminoAcid)chain.getAtomGroup(i + 1);
			//		System.out.println(rccto3d.Trans.getDihedral(a1, a2));

			//}
			//rccto3d.Trans.writePDB("datos/tst5.pdb", chain.getParent());

			/*
			rccto3d.Trans.makeBondTrans(chain, 0);
			rccto3d.Trans.makeBondTrans(chain, 1);
			rccto3d.Trans.makeBondTrans(chain, 2);
			rccto3d.Trans.makeBondTrans(chain, 3);

			rccto3d.Trans.writePDB("datos/tst1.pdb", chain.getParent());

			mueve(chain, angle, 0, true);
			mueve(chain, angle, 1, true);
			mueve(chain, angle, 2, true);
			mueve(chain, angle, 3, true);

			rccto3d.Trans.writePDB("datos/tst2.pdb", chain.getParent());
			*/

    } catch (Exception e) {
    	System.out.println("error");
      e.printStackTrace();
    }
  }

	public static void pega(Chain chain, double angle, int cual, boolean verbos) throws Exception{
		rccto3d.Trans.joinAmino(chain, cual);
		a1 = (AminoAcid)chain.getAtomGroup(cual);
		a2 = (AminoAcid)chain.getAtomGroup(cual + 1);
		double angle2 = rccto3d.Trans.getOmega(a1, a2);
		rccto3d.Trans.rotateOmega(chain, cual, -angle2*trans +angle*trans);
		angle2 = rccto3d.Trans.getOmega(a1, a2);
		if(verbos) System.out.println("ang " + angle2);
	}

	public static void mueve(Chain chain, double angle, int cual, boolean verbos) throws Exception{
		rccto3d.Trans.joinAmino(chain, cual);
		a1 = (AminoAcid)chain.getAtomGroup(cual);
		a2 = (AminoAcid)chain.getAtomGroup(cual + 1);

		double angle2 = rccto3d.Trans.getPhi(a1, a2);
		rccto3d.Trans.rotatePhi(chain, cual + 1, -angle2*trans +angle*trans);
		angle2 = rccto3d.Trans.getPhi(a1, a2);
		if(verbos) System.out.println("ang " + angle2);

		angle2 = rccto3d.Trans.getPsi(a1, a2);
		rccto3d.Trans.rotatePsi(chain, cual, angle2*trans -angle*trans);
		angle2 = rccto3d.Trans.getPsi(a1, a2);
		if(verbos) System.out.println("ang2 " + angle2);
	}

	public static void mueveAlfa(Chain chain, int cual, boolean verbos) throws Exception{
		rccto3d.Trans.joinAmino(chain, cual);
		a1 = (AminoAcid)chain.getAtomGroup(cual);
		a2 = (AminoAcid)chain.getAtomGroup(cual + 1);

		double angle2 = rccto3d.Trans.getPhi(a1, a2);
		rccto3d.Trans.rotatePhi(chain, cual + 1, -angle2*trans + A_PHI*trans);
		angle2 = rccto3d.Trans.getPhi(a1, a2);
		if(verbos) System.out.println("ang " + angle2);

		angle2 = rccto3d.Trans.getPsi(a1, a2);
		rccto3d.Trans.rotatePsi(chain, cual, angle2*trans -A_PSI*trans);
		angle2 = rccto3d.Trans.getPsi(a1, a2);
		if(verbos) System.out.println("ang2 " + angle2);
	}
	public static int[] calcRCC(String s1){
		int rcc[] = new int[26];
		try{
		String s2 = "A";
		Process p = Runtime.getRuntime().exec("python create_26dvRCC.py "+s1+" "+s2);
		BufferedReader in = new BufferedReader(new InputStreamReader(p.getInputStream()));
		String ret = in.readLine();
		String val = "";
		int j = 0;
		for(int i = 0; i < ret.length(); i++){
			char c = ret.charAt(i);
			if(c==','){
				rcc[j] = Integer.valueOf(val);
				val = "";
				j++;
			}else{
				val+=c;
			}
		}
		rcc[j] = Integer.valueOf(val);
		}catch(Exception e){
  		e.printStackTrace();
		}
		return rcc;
	}
}
