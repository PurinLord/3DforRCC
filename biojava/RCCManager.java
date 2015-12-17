package rccto3d;

import org.biojava.bio.structure.Structure;

import java.lang.Math;
import java.io.*;

public class RCCManager{

String fileDir;

public RCCManager(){
	this.fileDir = "./";
}

public RCCManager(String fileDir){
	this.fileDir = fileDir;
}

public int[] calcRCC(String s1, String tmpdir){
	int rcc[] = new int[26];
	try{
	String s2 = "A";
	Process p = Runtime.getRuntime().exec("python create_26dvRCC.py "+s1+" "+s2 +" "+tmpdir);
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

public int[] calcRCC(String s1){
	return calcRCC(s1, fileDir);
}

public int calcContactos(String s1, String tmpdir){
	int rcc[] = calcRCC(s1, tmpdir);
	int contactos = 0;
	int factor = 0;
	for(int i = 1; i <= rcc.length; i++){
		if(i <= 3){
			factor = 3;
		}else{
		if(i <= 5){
			factor = 4;
		}else{
		if(i <= 7){
			factor = 5;
		}else{
		if(i <= 11){
			factor = 6;
		}}}}
		contactos += rcc[i-1]*factor;
	}
	return contactos;
}

public int calcContactos(String s1){
	return calcContactos(s1, fileDir);
}

public double calcSimilarity(int[] rcc1, int[] rcc2){
	int sum = 0;
  int i=0;
	for(i=0; i<26; i++){
		sum += Math.pow(rcc1[i] - rcc2[i], 2.0);
	}
	return Math.sqrt(sum);
}// end calcSimilarity

public static void main(String args[]){
	RCCManager rccM = new RCCManager(args[1]);
	System.out.println(rccM.calcContactos(args[0]));
}
}
