package rccto3d;

import java.io.File;
import java.io.PrintWriter;
import java.lang.Math;
import java.lang.Double;

import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.Group;
import org.biojava.bio.structure.AminoAcid;
import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.AtomImpl;
import org.biojava.bio.structure.io.PDBFileReader;
import org.biojava.bio.structure.Calc;
import org.biojava.bio.structure.jama.Matrix;
import org.biojava.bio.structure.StructureException;

/** Makes a solid body rotatios of a spesific residue to change its dihedral angle
 */
public class Trans extends Calc {
  public static void main(String[] args) {
    Structure structure;
		Chain chain;
		AminoAcid a1;
		AminoAcid a2;
		double teta = 3.1416*(1.0/8);
		int res = 2;

    try {
			// SET PATH HERE!
			//structure = readPDB("1hiv.pdb");
			structure = readPDB(args[0]);
			chain = structure.getChain(0);

			a1 = (AminoAcid)chain.getAtomGroup(res);
			a2 = (AminoAcid)chain.getAtomGroup(res + 1);
			System.out.println(getDihedral(a1, a2));

			setPhi(chain, res, 0);
			setPsi(chain, res, 0);
			setOmega(chain, res, 0);

			System.out.println(getDihedral(a1, a2));

			//writePDB("modify0.pdb", structure);
			//printDihedral(chain, res);
			//rotatePhi(chain, res, teta);
			//printDihedral(chain, res);
			//writePDB("modify1.pdb", chain.getParent());

    } catch (Exception e) {
    	System.out.println("error");
      e.printStackTrace();
    }
  }

	private static final double BOND_DISTANCE = 1.32;
	private static final double OMEGA_TRANS = 180.0;//Math.PI;
	private static final double ANG_TO_RAD = Math.PI/180.0;
	private static final double X_BOND = 0.53689236886; //1.32 * Math.cos(180-114)
	private static final double Y_BOND = 1.20588000408; //1.32 * Math.sin(180-114)
	private static final double ANG_BOND = 123;

	/** Rotates residue resNumber of chain to change its dihedral angle Phi
	 */
	public static void rotatePsi(Chain chain, int resNumber, double angle){
		try{
			AminoAcid amino = (AminoAcid) chain.getAtomGroup(resNumber);
			//Stores original atom axis position
			Atom cAtom = (Atom) amino.getC().clone();
			double[] axis = getAxisPsi(chain, resNumber);
			Matrix rot = getRotMatrix(axis, angle);
			Calc.rotate(amino, rot);
			//Gets new atom axis position
			Atom rAtom = amino.getC();
			Atom tAtom = subtract(cAtom, rAtom);
			//Corrects
			Calc.shift(amino, tAtom);

			for(int i = resNumber - 1; i >= 0; i--){
				amino = (AminoAcid) chain.getAtomGroup(i);
				Calc.rotate(amino, rot);
				Calc.shift(amino, tAtom);
			}
		}catch (ClassCastException e){
    	e.printStackTrace();
			System.out.println("Chain has non-aminoacid elements");
		}catch (Exception e){
    	e.printStackTrace();
		}
	}

	/** Rotates residue resNumber of chain to change its dihedral angle Psi
	 */
	public static void rotatePhi(Chain chain, int resNumber, double angle){
		try{
			AminoAcid amino = (AminoAcid) chain.getAtomGroup(resNumber);
			Atom cAtom = (Atom) amino.getN().clone();
			double[] axis = getAxisPhi(chain, resNumber);
			Matrix rot = getRotMatrix(axis, angle);
			Calc.rotate(amino, rot);
			Atom rAtom = amino.getN();
			Atom tAtom = subtract(cAtom, rAtom);
			Calc.shift(amino, tAtom);

			for(int i = resNumber + 1; i < chain.getAtomLength(); i++){
				amino = (AminoAcid) chain.getAtomGroup(i);
				Calc.rotate(amino, rot);
				Calc.shift(amino, tAtom);
			}
		}catch (ClassCastException e){
    	e.printStackTrace();
			System.out.println("Chain has non-aminoacid elements");
		}catch (Exception e){
    	e.printStackTrace();
		}
	}

	public static void rotateOmega(Chain chain, int resNumber, double angle){
		try{
			AminoAcid amino = (AminoAcid) chain.getAtomGroup(resNumber);
			AminoAcid amino2 = (AminoAcid) chain.getAtomGroup(resNumber + 1);
			Atom cAtom = (Atom) amino2.getN().clone();
			double[] axis = getAxisOmega(chain, resNumber);
			Matrix rot = getRotMatrix(axis, angle);
			Calc.rotate(amino2, rot);
			Atom rAtom = amino2.getN();
			Atom tAtom = subtract(cAtom, rAtom);
			Calc.shift(amino2, tAtom);

		}catch (ClassCastException e){
    	e.printStackTrace();
			System.out.println("Chain has non-aminoacid elements");
		}catch (Exception e){
    	e.printStackTrace();
		}
	}

	public static void setPsi(Chain chain, int resNumber, double angle) throws Exception{
		AminoAcid a1 = (AminoAcid)chain.getAtomGroup(resNumber);
		AminoAcid a2 = (AminoAcid)chain.getAtomGroup(resNumber + 1);
		double angle2 = getPsi(a1, a2);
		rotatePsi(chain, resNumber, angle2*ANG_TO_RAD -angle*ANG_TO_RAD);
	}

	public static void setPhi(Chain chain, int resNumber, double angle) throws Exception{
		AminoAcid a1 = (AminoAcid)chain.getAtomGroup(resNumber);
		AminoAcid a2 = (AminoAcid)chain.getAtomGroup(resNumber + 1);
		double angle2 = getPhi(a1, a2);
		rotatePhi(chain, resNumber + 1, -angle2*ANG_TO_RAD +angle*ANG_TO_RAD);
	}

	public static void setOmega(Chain chain, int resNumber, double angle) throws Exception{
		AminoAcid a1 = (AminoAcid)chain.getAtomGroup(resNumber);
		AminoAcid a2 = (AminoAcid)chain.getAtomGroup(resNumber + 1);
		double angle2 = getOmega(a1, a2);
		rotateOmega(chain, resNumber, -angle2*ANG_TO_RAD +angle*ANG_TO_RAD);
	}

	public static void joinAmino(Chain chain, int resNumber) throws Exception{
		try{
			AminoAcid amino1 = (AminoAcid) chain.getAtomGroup(resNumber);
			AminoAcid amino2 = (AminoAcid) chain.getAtomGroup(resNumber + 1);
			Atom CA = amino1.getCA();
			Atom C = amino1.getC();
			Atom N2 = amino2.getN();
			Atom CA2 = amino2.getCA();

				Trans.writePDB("datos/tst0.pdb", chain.getParent());

			Atom e1 = subtract(C, CA);
			e1 = unitVector(e1);
				//System.out.println("e1 x " + e1.getX() + " y " + e1.getY() + " z " + e1.getZ());
			Atom e2 = makeAtom(-e1.getY(), e1.getX(), 0.0);
				//System.out.println("e2 x " + e2.getX() + " y " + e2.getY() + " z " + e2.getZ());
			Atom tAtom = add(add(C, scaleEquals(e1, X_BOND)), scaleEquals(e2, Y_BOND));
				//System.out.println("tAtom x " + tAtom.getX() + " y " + tAtom.getY() + " z " + tAtom.getZ());
			shift(amino2, tAtom);
				Trans.writePDB("datos/tst1.pdb", chain.getParent());

			Atom cAtom = (Atom) amino2.getN().clone();
				//System.out.println("cAtom x " + cAtom.getX() + " y " + cAtom.getY() + " z " + cAtom.getZ());
			double angle = angle(subtract(C, N2), subtract(CA2, N2));
				System.out.println("angle " + angle + " gama " + (ANG_BOND - angle));
			Atom axis = vectorProduct(subtract(C, N2), subtract(CA2, N2));
				System.out.println("axis " + axis.getZ());
			axis = unitVector(axis);
				System.out.println("axis " + axis.getZ());
			Matrix rot = getRotMatrix(new double[] {0,0,1}, -axis.getZ() * (ANG_BOND - angle)*ANG_TO_RAD);
			rotate(amino2, rot);
				Trans.writePDB("datos/tst2.pdb", chain.getParent());
				//Atom subtr = subtract(cAtom, N2);
				//System.out.println("subtr x " + subtr.getX() + " y " + subtr.getY() + " z " + subtr.getZ());
			shift(amino2, subtract(cAtom, N2));
				angle = angle(subtract(C, N2), subtract(CA2, N2));
				System.out.println("2angl " + angle);
				Trans.writePDB("datos/tst3.pdb", chain.getParent());

		}catch (ClassCastException e){
			e.printStackTrace();
			System.out.println("Chain has non-aminoacid elements");
		}catch (Exception e){
			e.printStackTrace();
		}
	}

	public static void makeBondTrans(Chain chain, int resNumber){
		try{
		joinAmino(chain, resNumber);
		//AminoAcid a1 = (AminoAcid)chain.getAtomGroup(resNumber);
		//AminoAcid a2 = (AminoAcid)chain.getAtomGroup(resNumber + 1);
		//rotateOmega(chain, resNumber, -getOmega(a1, a2)*ANG_TO_RAD + OMEGA_TRANS*ANG_TO_RAD);
		//setOmega(chain, resNumber, OMEGA_TRANS);
		}catch (Exception e){
			e.printStackTrace();
		}
	}

	public static double getOmega(AminoAcid a, AminoAcid b) throws StructureException {
		if ( ! Calc.isConnected(a,b)){
			throw new StructureException("can not calc Omega - AminoAcids are not connected!") ;
		}
		Atom a_CA = a.getCA();
		Atom a_C = a.getC();
		Atom b_N = b.getN();
		Atom b_CA = b.getCA();
		// C and N were checked in isConnected already
		if (b_CA==null) throw new StructureException("Can not calculate Omega, CA atom is missing");
		return torsionAngle(a_CA,a_C,b_N,b_CA);
	}

	/** Creates de axis of rotation to change dihedral angle Phi of residue ResNumber
	 */
	private static double[] getAxisPhi(Chain chain, int resNumber) throws Exception{
		AminoAcid amino = (AminoAcid) chain.getAtomGroup(resNumber);
		double[] axis;
		Atom a1 = amino.getN();
		Atom a2 = amino.getCA();
		axis = getDif(a1, a2);
		axis = normalize(axis);
		return axis;
	}

	/** Creates de axis of rotation to change dihedral angle Psi of residue ResNumber
	 */
	private static double[] getAxisPsi(Chain chain, int resNumber) throws Exception{
		AminoAcid amino = (AminoAcid) chain.getAtomGroup(resNumber);
		double[] axis;
		Atom a1 = amino.getCA();
		Atom a2 = amino.getC();
		axis = getDif(a1, a2);
		axis = normalize(axis);
		return axis;
	}

	private static double[] getAxisOmega(Chain chain, int resNumber) throws Exception{
		AminoAcid amino = (AminoAcid) chain.getAtomGroup(resNumber);
		AminoAcid amino2 = (AminoAcid) chain.getAtomGroup(resNumber + 1);
		double[] axis;
		Atom a1 = amino.getC();
		Atom a2 = amino2.getN();
		axis = getDif(a1, a2);
		axis = normalize(axis);
		return axis;
	}

	private static Atom getTransAtom(double[] vec){
		Atom atom = new AtomImpl();
		atom.setX(vec[0]);
		atom.setY(vec[1]);
		atom.setZ(vec[2]);
		return atom;
	}

	/** Returns the substraction of the coordinates of the atoms
	 */
	private static double[] getDif(Atom a1, Atom a2){
		double[] dif = new double[3];
		dif[0] = a1.getX() - a2.getX();
		dif[1] = a1.getY() - a2.getY();
		dif[2] = a1.getZ() - a2.getZ();

		return dif;
	}

	/** Constructs a rotation matrix of an angle teta along the axis vec
	 */
	public static Matrix getRotMatrix(double[] vec, double teta){
		double[][] rot = new double[3][3];
		double l = vec[0];
		double m = vec[1];
		double n = vec[2];
		rot[0][0] = l*l*(1.0 - Math.cos(teta)) + Math.cos(teta);
		rot[0][1] = m*l*(1.0 - Math.cos(teta)) - n*Math.sin(teta);
		rot[0][2] = n*l*(1.0 - Math.cos(teta)) + m*Math.sin(teta);
		rot[1][0] = l*m*(1.0 - Math.cos(teta)) + n*Math.sin(teta);
		rot[1][1] = m*m*(1.0 - Math.cos(teta)) + Math.cos(teta);
		rot[1][2] = n*m*(1.0 - Math.cos(teta)) - l*Math.sin(teta);
		rot[2][0] = l*n*(1.0 - Math.cos(teta)) - m*Math.sin(teta);
		rot[2][1] = m*n*(1.0 - Math.cos(teta)) + l*Math.sin(teta);
		rot[2][2] = n*n*(1.0 - Math.cos(teta)) + Math.cos(teta);
		return new Matrix(rot);
	}

	/** Calcula la lonfitud del vector
	 */
	private static double norm(double[] vec){
		double norm = Math.sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
		return norm;
	}

	/** Normalize vector
	 */
	public static double[] normalize(double[] vec){
		double norm = norm(vec);
		vec[0] = vec[0]/norm;
		vec[1] = vec[1]/norm;
		vec[2] = vec[2]/norm;
		return vec;
	}

	public static Atom makeAtom(double x, double y, double z){
		Atom a = new AtomImpl();
		a.setX(x);
		a.setY(y);
		a.setZ(z);
		return a;
	}

	public static String getDihedral(AminoAcid a1, AminoAcid a2){
		try{
		return getPhi(a1, a2) + " " + getPsi(a1, a2) + " " + getOmega(a1, a2);
		} catch (Exception e) {
		  e.printStackTrace();
			return "cant get angle";
		}
	}

	/** Creates a Structure objet from a PDB file
	 */
	public static Structure readPDB(String filename){
 		PDBFileReader pdbreader = new PDBFileReader();
		Structure structure = null;
    try{
    	structure = pdbreader.getStructure(filename);
    	//System.out.println(structure);
    } catch (Exception e) {
    	e.printStackTrace();
    }
    return structure;
	}

	/** Saves structure in to a pdb format in to file
	 */
	public static void writePDB(String filename, Structure structure){
		try{
		PrintWriter writer = new PrintWriter(filename, "UTF-8");
		writer.println(structure.toPDB().trim());
		writer.close();
		}catch(Exception e){
			e.printStackTrace();
		}
	}

	public static Atom scaleEquals(Atom a, double s) {
		double x = a.getX();
		double y = a.getY();
		double z = a.getZ();

		x *= s;
		y *= s;
		z *= s;

		//Atom b = new AtomImpl();
		a.setX(x);
		a.setY(y);
		a.setZ(z);

		return a;
	}

	/** Prints dhihedral angles Phi Psi between the two atoms
	 */
	private static void printPhiPsi(AminoAcid a1, AminoAcid a2){
		try{
			if(a1 != null && a2 != null){
				System.out.println(Calc.getPhi(a1, a2) + " " + Calc.getPsi(a1, a2));
			}
		}catch(Exception e){
			e.printStackTrace();
		}
	}

	/** Prints dhihedral near the given residue of chain
	 */
	public static void printDihedral(Chain chain, int resNumber){
		AminoAcid a1, a2, a3, a4, a5;
		a1 = a2 = a4 = a5 = null;
		try{
		if(resNumber >= 2){
			a1 = (AminoAcid)chain.getAtomGroup(resNumber - 2);
			if(resNumber >= 1){
				a2 = (AminoAcid)chain.getAtomGroup(resNumber - 1);
			}
		}
		a3 = (AminoAcid)chain.getAtomGroup(resNumber);
		if(chain.getAtomLength() - resNumber > 1){
			a4 = (AminoAcid)chain.getAtomGroup(resNumber + 1);
			if(chain.getAtomLength() - resNumber > 2){
				a5 = (AminoAcid)chain.getAtomGroup(resNumber + 2);
			}
		}
		System.out.println();	
		printPhiPsi(a1, a2);
		printPhiPsi(a2, a3);
		printPhiPsi(a3, a4);
		printPhiPsi(a4, a5);
    } catch (Exception e) {
      e.printStackTrace();
    }
	}

	public static double bondDist(Chain chain, int resNumber){
		try{
		AminoAcid amino1 = (AminoAcid) chain.getAtomGroup(resNumber);
		AminoAcid amino2 = (AminoAcid) chain.getAtomGroup(resNumber + 1);
		Atom atom1 = amino1.getC();
		Atom atom2 = amino2.getN();
		return getDistance(atom1, atom2);
		}catch(Exception e){
			e.printStackTrace();
			return 0.0;
		}
	}
}
