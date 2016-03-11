#include "protein2.h"

void Protein2::rdpdb(string fn){
	info_err = "";
	//FILE *fp = openfile(fn, "r");
	//if(fp == NULL) die("wrong file: %s\n", fn.c_str());
	char str[121], rn[8], an[5]; string rinfo0="", rinfo;
	double xt[3];
	name = basename((char*)fn.c_str());
  std::istringstream iss(fn);
  std::string line;    
  while (std::getline(iss,line)) {
		strcpy(str, line.c_str());
	//}
	//while(fgets(str,120,fp) != NULL){
		if(strstr(str, "END") == str) break;
		if(strstr(str, "TER") == str) break;
		if(strstr(str, "ATOM ") != str) continue;
		sscanf(str+13, "%3s", an); an[2] = '\0';
		if(strcmp(an, repAtom.c_str()) != 0) continue;
		sscanf(str+17, "%3s", rn); rn[3] = '\0';
		string line(str);
		rinfo = line.substr(17, 10);
		if(rinfo == rinfo0) continue;
//		str2dat(str+30, 3, xt);
		//<<xt coordenadas
		for(int m=0; m<3; m++) xt[m] = strtod(str+30+8*m, NULL);
		//std::cout << "xt " << xt[0]<<" "<<xt[1]<<" "<<xt[2] << "\n";
		x.push_back(Xvec(xt));
		//std::cout << "rn " << rn<< "\n";
		int id = aaDefine(rn, DEBUG);
		resid.push_back(id);
		sscanf(str+22, "%7s", rn); rn[7] = '\0';
		resnum.push_back(rn);
		rinfo0 = rinfo;
	}
	nres = resid.size();
	//fclose(fp);
	if(nres < 3) {
		info_err += " few residues";
		return;
	}
//
	calSS();
}
Protein2::Protein2(int nr, string seq1, double xp[][3]){
	info_err = "";
	name = "UNK";
	nres = nr;
	resid.resize(nr);
	for(int i=0; i<nr; i++) resid[i] = aaDefine1(seq1[i]);
	x.resize(nr);
	for(int i=0; i<nr; i++){
		x[i].setx(xp[i]);
	}
	if(nres < 3) {
		info_err += " few residues";
		return;
	}
	calSS();
}
void Protein2::calnneib(vector<int> &idx, double R0){
	assert(idx.size() == nres);
	if(R0 < 0){		// build idx for 3 distances
		calnneib0(idx); return;
	}
	for(int i=0; i<nres; i++){
		if(idx[i] == 1) continue;
		for(int j=0; j<nres; j++){
			if(idx[j] != 1) continue;
			double r2 = distance2(i, j);
			if(r2 < R0*R0){
				idx[i] = 2; break;
			}
		}
	}
	return;
}
void Protein2::calnneib0(vector<int> &idx){
	assert(idx.size() == nres);
	double r2_dim[] = {10., 12., 14.};
	for(int i=0; i<nres; i++){
		if(idx[i] == 1) continue;
		double rmin2 = 10000.;
		for(int j=0; j<nres; j++){
			if(idx[j] != 1) continue;
			double r2 = distance2(i, j);
			if(rmin2 < r2) continue;
			rmin2 = r2;
			if(rmin2 <= r2_dim[0]*r2_dim[0]) break;
		}
		for(int m=0; m<3; m++){
			if(rmin2 > r2_dim[m]*r2_dim[m]) continue;
			idx[i] = 2 + m; break;
		}
	}
}
void Protein2::calSS(){
	double Rdis[nres][3];
	const double Ra[] = {5.45, 5.18, 6.37}, Da = 2.1;
	const double Rb[] = {6.1, 10.4, 13}, Db = 1.42;
	if(nres <= 3) {
		fprintf(stderr, "%s: few residues %d", name.c_str(), nres);
		return;
	}
	ssec.resize(nres);
	for(int k=0; k<3; k++){
		for(int i=0; i<nres-k-2; i++) Rdis[i][k] = sqrt(distance2(i, i+k+2));
	}
	for(int i=0; i<nres; i++) ssec[i] = 0;
	for(int i=0; i<nres-5; i++){
		int iss = 0;
		if(fabs(Rdis[i][2]-Ra[2])<Da &&
			fabs(Rdis[i][1]-Ra[1])<Da && fabs(Rdis[i+1][1]-Ra[1])<Da &&
			fabs(Rdis[i][0]-Ra[0])<Da && fabs(Rdis[i+1][0]-Ra[0])<Da &&
			fabs(Rdis[i+2][0]-Ra[0])<Da) iss = 1;

		else if(fabs(Rdis[i][2]-Rb[2])<Db &&
			fabs(Rdis[i][1]-Rb[1])<Db && fabs(Rdis[i+1][1]-Rb[1])<Db &&
			fabs(Rdis[i][0]-Rb[0])<Db && fabs(Rdis[i+1][0]-Rb[0])<Db &&
			fabs(Rdis[i+2][0]-Rb[0])<Db) iss = 2;
		else continue;
//		for(int m=0; m<5; m++)
		ssec[i+2] = iss;
	}
// fragment
	vector<int> ifrag[2];
	int i1 = 0;
	for(int i=1; i<nres; i++){
		if(ssec[i]>0 && ssec[i-1]!=ssec[i]){
			ifrag[0].push_back(i); i1 = ssec[i];
		} else if(ssec[i] != i1){
			if(ifrag[1].size() < ifrag[0].size()) ifrag[1].push_back(i);
			i1 = ssec[i];
		}
	}
	if(ifrag[1].size() < ifrag[0].size()) ifrag[1].push_back(nres);
// remove short SSE & extend long SSE
	for(int i=0; i<ifrag[0].size(); i++){
		int i0=ifrag[0][i], i1=ifrag[1][i];
		if(i1-i0 >= 3) {
			if(i0 > 0) ssec[i0-1] = ssec[i0];
			if(i1 < nres-2) ssec[i1] = ssec[i1-1];
		} else {
			for(int j=i0; j<i1; j++) ssec[j] = 0;
		}
	}
// head & tail
	ssec[0] = ssec[1] = ssec[2];
	ssec[nres-1] = ssec[nres-2] = ssec[nres-3];
}
void Protein2::wrpdb(string fn){
	FILE *fp = openfile(fn, "w");
	wrpdb(fp);
	fclose(fp);
}
void Protein2::wrpdb(FILE *fp){
	char fmt1[] = "ATOM%7d  %-4s%3s %c%4d    %8.3f%8.3f%8.3f\n";
	for(int i=0; i<nres; i++){
		fprintf(fp, fmt1, i+1, "CA", rnam3_std[resid[i]], 'A',
				i, x[i][0], x[i][1], x[i][2]);
	}
	fprintf(fp, "TER\nEND\n");
}
double Protein2::distance2(int ia, int ib){
	return x[ia].distance2(x[ib]);
}
