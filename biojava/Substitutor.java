package rccto3d;

import java.util.Vector;
import java.util.Random;
import java.io.IOException;
import java.lang.Integer;
import java.util.ArrayList;
import java.util.Collections;

import org.biojava.bio.structure.io.PDBFileReader;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.AminoAcid;

import org.biojava.bio.structure.AminoAcid;

public class Substitutor{

Structure moldStructure = null;
Vector<Vector<Integer>> divition;
int divLength;
PDBfromFASTA pff = null;
Random rdm;

public Substitutor(Structure moldStructure){
	this.moldStructure = moldStructure;
	rdm =  new Random(System.currentTimeMillis());
}

public void randomInitialize(int divitionFactor){
	//Random rdm =  new Random(System.currentTimeMillis());
	int largo = moldStructure.getChain(0).getAtomLength()-1;
	int numSegment, maxSize, minSize, undefMax, undefMin;
	numSegment = rdm.nextInt((largo/divitionFactor)-4) + 4;
	Vector<Vector<Integer>> divition = new Vector<Vector<Integer>>(numSegment);
	Vector<Integer> segment;
	int start = 0;
	int size = 0;
	ArrayList<Integer> listStart = new ArrayList<Integer>(numSegment); 
	for(int i = 0; i < numSegment-1; i++){
		start = rdm.nextInt(largo-1)+1;
		listStart.add(start);
	}
	Collections.sort(listStart);
	start = 0;
	boolean par = rdm.nextBoolean();
	for(int i = 0; i < numSegment-1; i++){
		size = listStart.get(i) - start;
		if(par){
			segment = new Vector<Integer>(2);
			segment.add(start);
			segment.add(size);
			divition.add(segment);
		}
		par = !par;
		start = start + size;
	}
	size = largo - start;
	if(par){
		segment = new Vector<Integer>(2);
		segment.add(start);
		segment.add(size);
		divition.add(segment);
	}

	this.divLength = start + size;
	this.divition = divition;
}

public void randomInitialize(){
	int seed = rdm.nextInt(moldStructure.getChain(0).getAtomLength()/10);
	randomInitialize(seed);
}


@Deprecated
public void randomInitializeOld(){
	//Random rdm =  new Random(System.currentTimeMillis());
	int largo = moldStructure.getChain(0).getAtomLength();
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

public Vector<Vector<Integer>> getDivition(){
	return this.divition;
}

public void setDivition(Vector<Vector<Integer>> divition){
	this.divition = divition;
}

public void readDivition(String stringVector){
	Vector<Vector<Integer>> divition = new Vector<Vector<Integer>>();
	Vector<Integer> set = new Vector<Integer>(2);
	String[] split = stringVector.split(",");
	for(String a : split){
		if(set.size() >= 2){
			divition.add(set);
			set = new Vector<Integer>(2);
			set.add(Integer.parseInt(a.replace("[","").replace("]","").trim()));
		}else{
			set.add(Integer.parseInt(a.replace("[","").replace("]","").trim()));
		}
	}
	divition.add(set);
	this.divition = divition;
}

@Deprecated
public void createDivition(int numSegment, int maxSize, int minSize, int undefMax, int undefMin){
	//Vector<Integer>[] divition = (Vector<Integer>[]) new Vector<Integer>[numSegment];
	int largoStruc = this.moldStructure.getChain(0).getAtomLength();
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

/** Comienza la sutitución desde 0
 * TODO opciones para proteínas de longitudes distintas
 *
 */
public Structure fakeSubstitute(Structure subStructure){
	Chain chainFrom = moldStructure.getChain(0);
	Structure strucTo = (Structure)subStructure.clone();
	Chain chainTo = strucTo.getChain(0);
	AminoAcid a1;
	AminoAcid a2;
	double angle;
	Vector<Integer> segment;
	int index = 0;
	for(int i = 0; i < divition.size(); i ++){
		segment = divition.elementAt(i);
		for(int j = 0; j <= segment.elementAt(1); j++){
			try{
				index = segment.elementAt(0) + j;
				if(index+1 >= chainTo.getAtomLength()) break;
				a1 = (AminoAcid)chainFrom.getAtomGroup(index);
				a2 = (AminoAcid)chainFrom.getAtomGroup(index + 1);
				angle = Trans.getPhiFix(a1, a2);
				Trans.setPhi(chainTo, index, angle);
				angle = Trans.getPsiFix(a1, a2);
				Trans.setPsi(chainTo, index, angle);
				angle = Trans.getOmega(a1, a2);
				Trans.setOmega(chainTo, index, angle);
			}catch(Exception e){
  				e.printStackTrace();
			}
		}
	}
	return strucTo;
}

public Structure fakeSubstitute(Structure subStructure, int getFrom, int getTo, int setFrom){
	Chain chainFrom = moldStructure.getChain(0);
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
		//System.out.println("m m "+max+" "+min);
	if(max == min)return max;
	int size = rdm.nextInt(max-min);
	return size+min;
}

public static void main(String args[]){
	Substitutor sub = new Substitutor(Trans.readPDB(args[0]));
	sub.createDivition(Integer.parseInt(args[2]), 20, 10, 5, 2);
	Structure out = sub.fakeSubstitute(Trans.readPDB(args[1]));
	//sub.readDivition("[[4, 16], [26, 5], [33, 25], [67, 35], [110, 31]]");
	//Trans.writePDB("simDatos/out.pdb", out);
}

}
