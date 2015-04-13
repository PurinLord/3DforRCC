package rccto3d.managerot;

import java.io.*;
import java.util.*;
import java.lang.Math;

import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.Group;
import org.biojava.bio.structure.AminoAcid;
import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.AtomImpl;
import org.biojava.bio.structure.jama.Matrix;

import rccto3d.rotamers.*;
import rccto3d.*;

public class ArrangeAAPDB{

	private static final double RAD_TO_ANG = 180.0/Math.PI;
	private static final double ANG_TO_RAD = Math.PI/180.0;

	public static void main(String[] args) throws Exception{

		Structure struc;
		Chain chain;
		AminoAcid amino;
		Atom tAtom;
		Atom a1;
		Atom a2;
		Atom a3;
		double alpha;

		File directory = new File(args[0]);
		File[] fList = directory.listFiles();
		for (File file : fList){

			struc = Trans.readPDB(file.getPath());
			System.out.println(file.getName());
			struc = possition(struc);
			Trans.writePDB("aaPos/" + file.getName(), struc);
		}
	}

	public static Structure possition(Structure struc) throws Exception{

		Chain chain;
		AminoAcid amino;
		Atom tAtom;
		Atom a1;
		Atom a2;
		Atom a3;
		double alpha;

		chain = struc.getChain(0);
		amino = (AminoAcid)chain.getAtomGroup(0);

		tAtom = Trans.invert(amino.getN());

		Trans.shift(struc, tAtom);

		a1 = amino.getC();
		a2 = new AtomImpl();
		a2.setX(1.0);
		a2.setY(0.0);
		a2.setZ(0.0);
		a3 = Trans.vectorProduct(a1, a2);

		rotateTo(struc, a1, a2);

		a1 = (Atom) amino.getCA().clone();
		a1.setX(0.0);
		a2 = new AtomImpl();
		a2.setX(0.0);
		a2.setY(1.0);
		a2.setZ(0.0);
		a3 = Trans.vectorProduct(a1, a2);

		rotateTo(struc, a1, a2);

		return struc;
	}

	public static Structure rotateTo(Structure struc, Atom a1, Atom a2){
		Atom a3 = Trans.vectorProduct(a1, a2);
		double alpha = Trans.angle(a1, a2);
		double[] axis = new double[] {a3.getX(), a3.getY(), a3.getZ()};
		axis = Trans.normalize(axis);
		Matrix rot = Trans.getRotMatrix(axis, -alpha * ANG_TO_RAD);
		Trans.rotate(struc, rot);

		return struc;
	}
}
