package rccto3d;

import java.util.Vector;
import java.util.Random;
import java.io.IOException;
import java.lang.Integer;

import org.biojava.bio.structure.io.PDBFileReader;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.AminoAcid;

import org.biojava.bio.structure.AminoAcid;

public class Substitutor{

Structure subStructure = null;
Structure moldStructure = null;
Vector<Vector<Integer>> divition;
int divLength;
PDBfromFASTA pff = null;
Random rdm;

public Substitutor(Structure subStructure){
	this.subStructure = subStructure;
	rdm =  new Random(System.currentTimeMillis());
}

public void randomInitialize(){
	//Random rdm =  new Random(System.currentTimeMillis());
	int largo = subStructure.getChain(0).getAtomLength();
	int numSegment, maxSize, minSize, undefMax, undefMin;
	numSegment = rdm.nextInt(largo) + 2;
	maxSize = rdm.nextInt(largo) + 5;
	undefMax = rdm.nextInt(largo) + 1;
	if(undefMax > maxSize){int tmp = maxSize; maxSize = undefMax; undefMax = tmp;}//SWAP
	minSize = rdm.nextInt(maxSize);
	undefMin = rdm.nextInt(undefMax);
	int largoActual = numSegment*maxSize + numSegment*undefMax + undefMin;
	numSegment *= Math.sqrt(largo/(float)largoActual);
	maxSize *= Math.sqrt(largo/(float)largoActual);
	undefMax *= Math.sqrt(largo/(float)largoActual);
	undefMin *= Math.sqrt(largo/(float)largoActual);
	minSize *= Math.sqrt(largo/(float)largoActual);
	if(maxSize == 0){maxSize = 1;}
	if(undefMax == 0){undefMax = 1;}
	if(minSize == 0){minSize = 1;}
	if(undefMin == 0){undefMin = 1;}
	largoActual = numSegment*maxSize + numSegment*undefMax + undefMin;
	//System.out.println(numSegment+"\t"+maxSize+"\t"+minSize+"\t"+undefMax+"\t"+undefMin+"\t"+largoActual);

	this.createDivition(numSegment, maxSize, minSize, undefMax, undefMin);
}

public void setDivition(Vector<Vector<Integer>> divition){
	this.divition = divition;
}

public Vector<Vector<Integer>> getDivition(){
	return this.divition;
}

public void createDivition(int numSegment, int maxSize, int minSize, int undefMax, int undefMin){
	//Vector<Integer>[] divition = (Vector<Integer>[]) new Vector<Integer>[numSegment];
	int largoStruc = this.subStructure.getChain(0).getAtomLength();
	//int largoDiv = (numSegment * maxSize) + (numSegment * undefMax);
	//if(largoStruc < largoDiv){
	//	throw new IllegalArgumentException();
	//}
	Vector<Vector<Integer>> divition = new Vector<Vector<Integer>>(numSegment);
	Vector<Integer> segment;
	int size = 0;
	int start = 0;
	for(int i = 0; i < numSegment; i++){
		start += size + generateRandom(undefMax, undefMin);
		size = generateRandom(maxSize, minSize);
		segment = new Vector<Integer>(2);
		segment.add(start);
		segment.add(size);
		divition.add(segment);
	}
	this.divLength = start + size;
	this.divition = divition;
}

public Structure fakeSubstitute(Structure moldStructure){
	Chain chainFrom = moldStructure.getChain(0);
	Chain chainTo = subStructure.getChain(0);
	AminoAcid a1;
	AminoAcid a2;
	double angle;
	Vector<Integer> segment;
	int index = 0;
	for(int i = 0; i < divition.size(); i ++){
		segment = divition.elementAt(i);
		for(int j = 0; j < segment.elementAt(1); j++){
			try{
				index = segment.elementAt(0) + j;
				//System.out.println(index);
				a1 = (AminoAcid)chainFrom.getAtomGroup(i);
				a2 = (AminoAcid)chainFrom.getAtomGroup(i + 1);
				angle = Trans.getPhiFix(a1, a2);
				Trans.setPhi(chainTo, i, angle);
				angle = Trans.getPsiFix(a1, a2);
				Trans.setPsi(chainTo, i, angle);
				angle = Trans.getOmega(a1, a2);
				Trans.setOmega(chainTo, i, angle);
			}catch(Exception e){
  			e.printStackTrace();
			}
		}
	}
	return subStructure;
}

public Structure fakeSubstitute(Structure fitStruc, int getFrom, int getTo, int setFrom){
	Chain chainFrom = fitStruc.getChain(0);
	Chain chainTo = subStructure.getChain(0);
	double angle;
	for(int i = 0; i < getTo - getFrom-1; i ++){
		try{
			angle = Trans.getPhiFix((AminoAcid)chainFrom.getAtomGroup(i+getFrom), (AminoAcid)chainFrom.getAtomGroup(i+getFrom+1));
			Trans.setPhi(chainTo, i+setFrom, angle);
			angle = Trans.getPsiFix((AminoAcid)chainFrom.getAtomGroup(i+getFrom), (AminoAcid)chainFrom.getAtomGroup(i+getFrom+1));
			Trans.setPsi(chainTo, i+setFrom, angle);
		}catch(Exception e){
  		e.printStackTrace();
		}
	}
	return subStructure;
}

public int generateRandom(int max, int min){
	if(max == min)return max;
	int size = rdm.nextInt(max-min);
	return size+min;
}

public static void main(String args[]){
	Substitutor sub = new Substitutor(Trans.readPDB(args[0]));
	sub.createDivition(Integer.parseInt(args[2]), 20, 10, 5, 2);
	Structure out = sub.fakeSubstitute(Trans.readPDB(args[1]));
	Trans.writePDB("simDatos/out.pdb", out);
}

}
