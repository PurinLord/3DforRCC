#include "protein2.h"

const double GDT_r2[] = {1, 2*2, 4*4, 8*8};
inline double getd0_TM(double L0){
	if(L0 > 15) return max(.5, 1.24 * pow(L0-15., 1./3) - 1.8);
	return 0.5;
}
const double Eff0 = 1.3;	// to include those pairs between 8~9A
double Salign::Align_DP(){
    double dmax = -100.;
    
    if (nonSeqAlign)
        Run_NonSeq();
    else {
        double dat[3];
        for(int j=0; j<nB+1; j++) {Smat[0][j] = 0; Idir[0][j] = -1;}
        for(int i=0; i<nA+1; i++) {Smat[i][0] = 0; Idir[i][0] = -1;}
        Idir[0][0] = 0;
        double dmax = -100.;
        int imax=-1, jmax=-1;
        for(int i=1; i<nA+1; i++)
        for(int j=1; j<nB+1; j++){
            double dt = Rmat[i-1][j-1];
            dat[0] = Smat[i-1][j-1] + dt;
            dat[1] = Smat[i][j-1];
            dat[2] = Smat[i-1][j];
            if(Idir[i][j-1] != 1) dat[1] -= gap0;		// simple DP
            else dat[1] -= gap1;
            if(Idir[i-1][j] != 2) dat[2] -= gap0;
            else dat[2] -= gap1;
    //
            int it = 0;
            if(dat[it] < dat[1]) it = 1;
            if(dat[it] < dat[2]) it = 2;
            Smat[i][j] = dat[it]; Idir[i][j] = it;
            if(dat[it] < 0) { Smat[i][j] = 0; Idir[i][j] = -1; }
            if(dmax < Smat[i][j]) {dmax=Smat[i][j]; imax=i; jmax=j;}
        }
    // trace
        int imax0 = imax, jmax0 =jmax;
        if(imax<0 || jmax<0){
            printf("%d %d %f\n", imax, jmax, dmax);
            die("not used smat: %s %s", pa->name.c_str(), pb->name.c_str());
        }
        for(int m=0; m<2; m++) ialign[m].clear();
        int i0, j0, k0;
        while(imax>=0 || jmax>=0){
            if(imax==0 && jmax==0) break;
            k0 = Idir[imax][jmax];
            i0 = imax - 1; j0 = jmax - 1;
            if(k0 == 0) {imax --; jmax --;}
            else if(k0 == 1) {jmax --; i0 = -1;}
            else if(k0 == 2) {imax --; j0 = -1;}
            else if(k0 == -1) break;
            else die("not known k0: %d\n", k0);
            ialign[0].push_back(i0);
            ialign[1].push_back(j0);
        }
        for(int m=0; m<2; m++) reverse(ialign[m].begin(), ialign[m].end());
        return dmax;
    //
    // add the tail
        i0 = imax0; j0 = jmax0;
        while(i0<nA || j0<nB){
            if(i0 >= nA) ialign[0].push_back(-1);
            else ialign[0].push_back(i0);
            if(j0 >= nB) ialign[1].push_back(-1);
            else ialign[1].push_back(j0);
            i0 ++; j0 ++;
        }
        return dmax;
    }
    return dmax;
}
int Salign::redo_fit(int flag, double rcut2=64){
	using namespace Salign_PARAMS;
// flag: 0: all; 1: w=(1+r2/d2); 2:
	double xc1[3], xc2[3];
	if(flag == 1 && rcut2 < 0.1) rcut2 = D0_2;
	bool bSP = (score_type == iSP);
	bool bGDT = (score_type == iGDT);
	int nfit = 0;
	vector<int> list2[2];
	for(int m=0; m<2; m++) list2[m].clear();
	for(int i=0; i<ialign[0].size(); i++){
		int i1=ialign[0][i], i2=ialign[1][i];
		if(i1<0 || i2<0) continue;
		if(flag == 0) {
			wfit[nfit] = 1;
		}else {
			double r2 = xn1[i1].distance2( (*xbp)[i2] );
			if(flag == 1) {
				if((bSP ||bGDT) && r2>4*D0_2*Eff0) continue;
				double dt = r2 / rcut2;
				wfit[nfit] = 1 / (1 + dt);
//				wfit[nfit] = 1 / pow(1 + dt,2); // test on 03/21/12, <0.1% diff
			} else if(flag == 2){
				if(r2 > rcut2) continue;
				wfit[nfit] = 1.;
			}
		}
// to protect xn1 from being covered by small nfit
		list2[0].push_back(i1);
		list2[1].push_back(i2);
		nfit ++;
	}
	if(nfit < 3){
		if(DEBUG > 0) fprintf(stderr, "too small alignment size: %lu %d\n", ialign[0].size(), nfit);
		return 1;
		rotmol(nA, u, (*xap), *xap);
		pa->wrpdb("1.pdb"); pb->wrpdb("2.pdb");
		exit(0);
	}
//
	assert(nfit == list2[0].size());
	for(int i=0; i<nfit; i++){
		int i1 = list2[0].at(i), i2 = list2[1][i];
		for(int m=0; m<3; m++){
			xn1[i][m] = (*xap)[i1][m];
			xn2[i][m] = (*xbp)[i2][m];
		}
	}
	double ds = 0.;
	for(int i=0; i<nfit; i++) ds += wfit[i];
	if(ds > 1.0e-8){
		for(int i=0; i<nfit; i++) wfit[i] /= ds;
	} else {
		for(int i=0; i<nfit; i++) wfit[i] = 1./nfit;
	}
	translate2(nfit, wfit, xn1, xn1, xc1);
	translate2(nfit, wfit, xn2, xn2, xc2);
	quatfit(nfit, wfit, xn2, xn1, u);
	for(int m=0; m<3; m++) u[3][m] = xc2[m] - dot_product(u[m], xc1);
	rotmol(nA, u, (*xap), xn1);
	return 0;
}
int Salign::calrms(double &rms, double cut2=-1.){
	int n1 = 0; rms = 0.;
	for(int i=0; i<ialign[0].size(); i++){
		int i1=ialign[0][i], i2=ialign[1][i];
		if(i1<0 || i2<0) continue;
		double r2 = xn1[i1].distance2( (*xbp)[i2] );
		if(cut2 > 0 && r2 > cut2) continue;
		rms += r2;
		n1 ++;
	}
	if(n1 > 0) rms = sqrt(rms / n1);
	return n1;
}
// fill Rmat for DP by using "xn1, xbp, D0_2"
void Salign::calRmat(){
	for(int i=0; i<nA; i++)
	for(int j=0; j<nB; j++){
		double r2 = xn1[i].distance2((*xbp)[j]);
		if(score_type == iSP){
			Rmat[i][j] = 0.;
			if(r2 < D0_2*4*Eff0) Rmat[i][j] = 1. / (1. + r2 / D0_2) - 1/(1+4*Eff0);		// slightly extended to include those 8~9A
		} else if(score_type == iGDT) {
			Rmat[i][j] = 0.;
			if(r2 > GDT_r2[3]){
				if(r2 < GDT_r2[3]*Eff0) Rmat[i][j] = 1.e-3;
				continue;
			}
			for(int m=3; m>=0; m--){
				if(r2 > GDT_r2[m]) break;
				Rmat[i][j] = 1. - m/4.;
			}
		} else {
			Rmat[i][j] = 1. / (1. + r2 / D0_2);
		}
	}
}
double Salign::optimize_score(int i0=0, int ie=6){
	double dmax1 = 0., d1;
	double r2 = powi(D0+2., 2);
//	int types[] = {0,2,2,1}; int rs[] = {0,r2+1,r2,0};
	int types[] = {0,1,2,2,1}; int rs[] = {0,0,r2+1,r2,0};
	for(int m=i0; m<ie; m++){
		int ierr = 0;
		if(m >= 4) ierr = redo_fit(1,0);
		else ierr = redo_fit(types[m], rs[m]);
		if(ierr > 0 && m == 1){			// rare event in SP/GDT
			double dt = r2;
			for(int k=0; k<20; k++){
				dt *= 1.5;
				ierr = redo_fit(2, dt);
				if(ierr == 0) break;
			}
			if(ierr !=0) return dmax1;
		}
		if(ierr > 0) continue;
		if(ie>1 && m==0) continue;		// 2nd is always better
		d1 = calRscore();
		if(m>=4 && d1<dmax1) break;	// mature break after 4 turns
		if(dmax1 < d1) dmax1 = d1;
	}
	return dmax1;
}
double Salign::optimize_align(){
	gap0 = gap1 = 0.;
	double dmax1=-100.;
	double gap0_dim[] = {0.5, 0};
	for(int k=0; k<2; k++){
		rotmol(nA, u_sv, (*xap), xn1);
		gap0 = gap0_dim[k];
		for(int i=0; i<9; i++){
			calRmat();
			Align_DP();
			double d1 = optimize_score(1,20);
			if(d1-dmax1 < 0.01*dmax1) break;
			if(dmax1 < d1) dmax1 = d1;
		}
	}
	return dmax1;
}
double Salign::run_Align1(){
	int ierr = 0;
	gap0 = gap1 = 0.;
	calRmat();
	Align_DP();
	return optimize_score(1);
}
void Salign::Run_all(){
	if(bscoreOnly >= 1){
		scoreOnly(); return;
	}
//
	Rscore_max = -100.;
	if(bcheck > 0) {
		Run_custom();
		if(bcheck == 1) return;
	}
	Run1();
	if(Rscore_max >= 1.) return;
	Run2();
	Run3();
	Run4();
//	score_final();
}
void Salign::Run_custom(){
	FILE *fp = openfile(fali);
	char str[201], ss[3][20];
	while(fgets(str, 200, fp) != NULL){
		str2dat(str, 2, ss);
		ialign[0].push_back(atoi(ss[0]));
		ialign[1].push_back(atoi(ss[1]));
	}
	optimize_score();
	optimize_align();
}
// init seed from gapless threading
void Salign::Run1(){
	Rscore_sv = -100.;
	int nh = nA / 2;
	double scores[nB]; bzero(scores, sizeof(scores));
	for(int k=-nh; k<nB-nh; k++){
		for(int m=0; m<2; m++) ialign[m].clear();
		for(int j=0; j<nB; j++){
			if(j-k<0) continue;
			if(j-k>=nA) break;
			ialign[0].push_back(j-k);
			ialign[1].push_back(j);
		}
		if(ialign[0].size() < 3) continue;
		scores[k+nh] = optimize_score(0,2);
	}
	if(Rscore_max >= 1.) return;
//
	double dmax1 = -100.;
	for(int m=0; m<10; m++){
		int imax = maxloc(nB, scores);
		if(scores[imax] <= 0) break;
		if(m == 0) dmax1 = scores[imax];
		if(scores[imax] <= 0.8*dmax1) break;
		for(int m=0; m<2; m++) ialign[m].clear();
		int k = imax - nh;
		for(int j=0; j<nB; j++){
			if(j-k<0) continue;
			if(j-k>=nA) break;
			ialign[0].push_back(j-k);
			ialign[1].push_back(j);
		}
		if(ialign[0].size() < 3) continue;
		redo_fit(0);
		double d1 = run_Align1();
		for(int i=-5; i<=5; i++){		// del around to save time
			int i1 = imax + i; 
			if(i1>=0 && i1<nB) scores[i1] = -100.;
		}
	}
	optimize_align();
}
// init seed from secondary structures
void Salign::Run2(){
	Rscore_sv = -100.;
	gap0 = 3; gap1 = 0.;
	for(int i=0; i<nA; i++)
	for(int j=0; j<nB; j++){
		int i1 = pa->ssec[i], j1 = pb->ssec[j];
		double dt = 0.;
		if(i1 == j1) dt = 1.;
		else dt = 0.;
		Rmat[i][j] = dt;
	}
//
	double gap0_dim[] = {5,1,0};
	for(int m=0; m<3; m++){
		gap0 = gap0_dim[m];
		Align_DP();
		optimize_score();
	}
	optimize_align();
}
// init seed by fragment
void Salign::Run3(){
	Rscore_sv = -100.;
	int nmin = min(nA, nB);
	int nfrag = max(20, min(50, nmin/5));
	if(fragsize >= 1) nfrag = max(fragsize, 5);
	if(nfrag > nmin) nfrag = nmin;
	double rms1;
	int n0 = nfrag, n1 = nfrag;
	double rms0 = pow(nfrag, 1./3) * 1.2;
	for(int i=0; i<nA-nfrag+1; i++)
	for(int j=0; j<nB-nfrag+1; j+=n1){
		for(int m=0; m<2; m++) ialign[m].clear();
		for(int m=0; m<nfrag; m++){
			ialign[0].push_back(i+m);
			ialign[1].push_back(j+m);
		}
		redo_fit(0);
		calrms(rms1);
		if(rms1 > rms0) continue;
		for(int m=0; m<3; m++) redo_fit(1, 1.);
		double d1 = run_Align1();
	}
	optimize_align();
}
// init seed by DP on previous alignment + secondary structure
void Salign::Run4(){
	Rscore_sv = -100.;
	rotmol(nA, u_max, (*xap), xn1);
	calRmat();
	for(int i=0; i<nA; i++)
	for(int j=0; j<nB; j++){
		int i1 = pa->ssec[i], j1 = pb->ssec[j];
		if(i1 == j1) Rmat[i][j] += 0.5;
	}
	gap0 = gap1 = 0.;
	double gap0_dim[] = {5, 1., 0};
	for(int m=0; m<3; m++){
		gap0 = gap0_dim[m];
		Align_DP();
		optimize_score();
	}
	optimize_align();
}
// calculate score based on values of xn1/xbp; ialign & D0_2
double Salign::calRscore(){
	Rscore = 0.;
	for(int i=0; i<ialign[0].size(); i++){
		int i1 = ialign[0][i], i2 = ialign[1][i];
		if(i1 < 0 || i2 < 0) continue;
		double r2 = xn1[i1].distance2((*xbp)[i2]);
		if(score_type == iSP) {
			if(r2 > D0_2*4) continue;
			Rscore += 1.25 * (1./ (1. + r2/D0_2) - 0.2);
		} else if(score_type == iGDT) {
			double dt = 0.;
			for(int m=3; m>=0; m--){
				if(r2 > GDT_r2[m]) break;
				dt = 1 - m/4.;
			}
			Rscore += dt;
		} else {
			Rscore += 1. / (1. + r2/D0_2);
		}
	}
	Rscore /= Lmin;
	save_align();
	return Rscore;
}
void Salign::calscores_all(){
	double nali, rms, GDT, SP0, LG0, nid1;
	SP0 = LG0 = GDT = nali = rms = nid1 = 0.;
	double D2 = D00*D00;
// TMs
	double TMs[3], lTM[3], DTM2[3];
	bzero(TMs, sizeof(TMs));
	lTM[0] = (nA + nB) * 0.5;
	lTM[1] = min(nA, nB);
	lTM[2] = max(nA, nB);
	for(int m=0; m<3; m++) {
		double dt = getd0_TM(lTM[m]);
		DTM2[m] =  dt*dt;
	}
//
	double r2 = 0;
	int i1, i2 = 0;
		//printf("&%d %d\n",ialign[0].size(),ialign[1].size());
	int length = min(nA, nB);
	for(int i=0; i<length; i++){
		i1 = ialign[0][i];
		i2 = ialign[1][i];
			//printf("i- %d %d\n",i1,i2);
		//if(i1 < 0 || i2 < 0) continue;
		r2 = xn1[i1].distance2((*xbp)[i2]);
		for(int m=0; m<3; m++) TMs[m] += 1. / (1. + r2/DTM2[m]);
		LG0 += 1. / (1. + r2/D2);
			//printf("r2>%.2f\n",r2);
		for(int m=0; m<4; m++){
				//printf("<%.2f\n",GDT_r2[m]);
			if(r2 > GDT_r2[m]) continue;
			GDT += 1 - m/4.; break;
		}
		if(r2 < D2*4) SP0 += 1.25 * (1/ (1. + r2/D2) - 0.2);
		//if(r2 < 64.) { nali ++; rms += r2; }
		nali ++; rms += r2;
		if(pa->resid[i1] == pb->resid[i2]) nid1 ++;
			//printf("/ %d %.2f %.2f\n",i ,rms, nali);
	}
	scores_all[ieLA] = nali;
	nali = max(1., nali);
	rms = sqrt(rms / nali);
	scores_all[ieRMS] = rms;
	scores_all[ieSP] = SP0 / Lmin;
	if(score_type == iSP){
		int num[2][4]; calLali2(&num[0][0], Denv);
		assert(num[0][0] == num[1][0]);
		scores_all[ieLE] = (num[0][1] + num[1][1])*0.5 + num[0][0];
	}
	for(int m=0; m<3; m++){
		TMs[m] /= lTM[m];
	}
	scores_all[ieTMa] = TMs[0];
	scores_all[ieTMb] = TMs[1];
	scores_all[ieTMc] = TMs[2];
	scores_all[ieGDT] = GDT / Lmin;
	scores_all[ieLG] = LG0 / Lmin;
	scores_all[ieSEQ] = 100.*nid1 / Lmin;
	return;
}
bool Salign::save_align(){
	if(Rscore_sv > Rscore) return -1;
	Rscore_sv = Rscore;
	for(int m=0; m<12; m++) u_sv[0][m] = u[0][m];
	for(int m=0; m<2; m++) ialign_sv[m] = ialign[m];
//
	if(Rscore_max > Rscore) return 0;
	Rscore_max = Rscore;
	for(int m=0; m<12; m++) u_max[0][m] = u[0][m];
	for(int m=0; m<2; m++) ialign_max[m] = ialign[m];
	return 1;
}
void Salign::print_max(FILE *fp){
	string sinfo;
	print_max(sinfo);
	fprintf(fp, "%s", sinfo.c_str());
	fflush(fp);
}
void Salign::print_max(string &sinfo){
	for(int m=0; m<12; m++) u[0][m] = u_max[0][m];
	for(int m=0; m<2; m++) ialign[m] = ialign_max[m];
	score_final();
	prtali(sinfo);
}
void Salign::pro_swap(int bmat=0){
	breverse = (! breverse);
	Protein2 *pt = pa; pa = pb; pb = pt;
	nA = pa->nres; nB = pb->nres;
	xap = &pa->x; xbp = &pb->x;
	if(! bmat) return;
// swap the matrix & ialign
	for(int i=0; i<3; i++)
	for(int j=i+1; j<3; j++){
		double dt = u[i][j];
		u[i][j] = u[j][i]; u[j][i] = dt;
	}
	double dat[3];
	for(int m=0; m<3; m++) dat[m] = -dot_product(u[m], u[3]);
	for(int m=0; m<3; m++) u[3][m] = dat[m];
	vector<int> vt = ialign[0];
	ialign[0] = ialign[1]; ialign[1] = vt;
}
void Salign::init(Protein2 *a, Protein2 *b, int brev){		// nA < nB
	Rscore_max = Rscore_sv = -100;
	pa = a; pb = b; breverse = 0;
	if(brev && pa->nres > pb->nres){
		pb = a; pa = b; breverse = 1;
	}
	nA = pa->nres; nB = pb->nres;
	xap = &pa->x; xbp = &pb->x;
	if(nA < 1 || nB < 1) die("one empty file\n");
	Lmin = min(nA, nB);
	D0 = D00;
	if(score_type == iTM) D0 = getd0_TM(double(Lmin));
	D0_2 = D0 * D0;
//
	if(nA>ma || nB>mb) alloc_x(); if(DEBUG > 0) printf("Nres: %s %d -- %s %d; cutoff: %.1f\n", 
		pa->name.c_str(), nA, pb->name.c_str(), nB, D0);
}
// swap back if in reverse status, and compute all scores for Rscore_max
void Salign::score_final(){
	if(breverse) pro_swap(1);
	rotmol(nA, u, (*xap), xn1);
	if (nonSeqAlign) Refine_NonSeq();
	calscores_all();
}
void Salign::calLali2(int *num, double R0){
	vector<int> idx[2];
	idx[0].resize(nA); idx[1].resize(nB);
	for(int i=0; i<nA; i++) idx[0][i] = -1;
	for(int i=0; i<nB; i++) idx[1][i] = -1;
//
	for(int i=0; i<ialign[0].size(); i++){
		int i1 = ialign[0][i], i2 = ialign[1][i];
		if(i1 < 0 || i2 < 0) continue;
		double r2 = xn1[i1].distance2((*xbp)[i2]);
		if(r2 > 64.) continue;
		idx[0][i1] = idx[1][i2] = 1;		// core region
	}
	pa -> calnneib(idx[0], R0);
	pb -> calnneib(idx[1], R0);
	for(int m=0; m<8; m++) num[m] = 0;
	for(int k=0; k<2; k++){
		for(int i=0; i<idx[k].size(); i++){
			int it = idx[k][i] - 1;
			if(it < 0 || it > 3) continue;
			num[k*4 + it] ++;
		}
	}
}
void Salign::delete_x(){
	if(ma<1 && mb<1) return;
}
void Salign::alloc_x(){
	delete_x();
	ma = nA; mb = nB;
	int nmax = max(ma, mb);
	wfit.resize(nmax, 1.);
	xn1.resize(nmax); xn2.resize(nmax);
	build_2D_array(Smat, ma+1, mb+1);
	build_2D_array(Rmat, ma+1, mb+1);
	build_2D_array(Idir, ma+1, mb+1);
}
void Salign::prtali(FILE *fp){
	string sinfo;
	prtali(sinfo);
	fprintf(fp, "%s", sinfo.c_str());
}
void Salign::calSPscores(double *sp){
	double sp0 = scores_all[ieSP]*Lmin;
	double L1 = max(scores_all[ieLE], 3.);
	sp[0] = sp0 / (pow(L1, 1.-Alpha) * 3.75);
	sp[1] = sp0 / (pow((nA+nB)*0.5, 1.-Alpha) * 3.75);
	sp[2] = sp0 / (pow(Lmin, 1.-Alpha) * 3.75);
	sp[3] = 1. / (1. + exp(-(sp[0] - 0.523) / 0.044));
}
void Salign::prtali(string &sinfo){
	if(rmsOnly){
		char str[20];
		//sprintf(str,"%.2f", scores_all[ieRMS]);
		printf("%.2f", scores_all[ieRMS]);
	}else{
	sinfo = "";
	char str[max(nA+nB+100,2001)];
//
	double SPs[4];		// SPe, SPa, SPb
	if(score_type == iSP) calSPscores(SPs);
//
	if(iprint >= 99) {
		sprintf(str, "################## Alignment Report ##################\n"); sinfo += str;
		sprintf(str, "%s %s: Length= %d %d\n", pa->name.c_str(), pb->name.c_str(), nA, nB); sinfo += str;
		if(score_type == iSP){
			double p0 = 100 * SPs[3];
			sprintf(str, "Pfold= %.1f %%; SPe/SPa/SPb= %.3f %.3f %.3f ;Effective_Length: %d\n", p0, SPs[0], SPs[1], SPs[2], int(scores_all[ieLE])); sinfo += str;
		}
		sprintf(str, "RMSD/Nali= %.2f / %d ;GDT= %.3f ;TMscore(a,b,c)= %.3f %.3f %.3f SEQID= %.1f%%\n", scores_all[ieRMS], int(scores_all[ieLA]+0.1), scores_all[ieGDT], scores_all[ieTMa], scores_all[ieTMb], scores_all[ieTMc], scores_all[ieSEQ]); sinfo += str;
		sprintf(str, "\nRotation Matrix:\n"); sinfo += str;
	} else {
		double e1 = scores_all[score_type];
		if(score_type == iSP) e1 = SPs[0];
		if(e1 < cutoff) return;
//
		sprintf(str, "%s %s ", pa->name.c_str(), pb->name.c_str()); sinfo+=str;
		if(score_type == iSP){
			sprintf(str, "%.3f %.3f %.3f %.3f %d %.1f\n", SPs[0], SPs[1], SPs[2], SPs[3], int(scores_all[ieLE]), scores_all[ieSEQ]);
		} else if(score_type == iTM){
			sprintf(str, "%.3f %.3f %.3f\n", scores_all[ieTMa], scores_all[ieTMb], scores_all[ieTMc]);
		} else {
			sprintf(str, "%.3f\n", scores_all[score_type]);
		}
		sinfo += str;
	}
	if(iprint <= 1) return;
//
// print matrix
	for(int m=0; m<3; m++){
		sprintf(str, "%10.5f %8.5f %8.5f %8.5f\n", u[3][m], u[m][0], u[m][1], u[m][2]); sinfo += str;
	}
	if(iprint == 2) return;
//
	string sa1, sa2, sdis; sa1 = sa2 = sdis = "";
	int SOres = 0; double SO; // Structure Overlap: defined as the percentage of aligned representative
	                          // atoms within 3.5Å of corresponding superimposed atoms.
	for(int i=0; i<ialign[0].size(); i++){
		int i1=ialign[0][i], i2=ialign[1][i];
		if(i1 < 0) sa1 += '-';
		else sa1 += rnam1_std[pa->resid[i1]];
		if(i2 < 0) sa2 += '-';
		else sa2 += rnam1_std[pb->resid[i2]];
		char c1 = ' ';
		if(i1>=0 && i2>=0) {
			double r2 = xn1[i1].distance2( (*xbp)[i2] );
			if(r2 <= 12.25) SOres++;
			if(r2 <= 16) c1 = ':';
			else if(r2 < 64) c1 = '.';
		}
		sdis += c1;
	}
	SO = (SOres / (double)min(nA,nB))*100;
//
	int is1=0, is2=0, ie1=0, ie2=0;
	for(int i=0; i<ialign[0].size(); i++){
		if(ialign[0][i] >= 0) {is1=ialign[0][i]; break;}
	}
	for(int i=ialign[0].size()-1; i>=0; i--){
		if(ialign[0][i] >= 0) {ie1=ialign[0][i]; break;}
	}
	for(int i=0; i<ialign[0].size(); i++){
		if(ialign[1][i] >= 0) {is2=ialign[1][i]; break;}
	}
	for(int i=ialign[1].size()-1; i>=0; i--){
		if(ialign[1][i] >= 0) {ie2=ialign[1][i]; break;}
	}
	if(iprint >= 99){
		sprintf(str, "\nAlignment: %s %s\t\tStructure Overlap: %.2f %%\n", pa->name.c_str(), pb->name.c_str(), SO); sinfo += str;
		sprintf(str, "(':' denotes the residue pairs of distance <= 4A, and '.' denotes <=8A)\n"); sinfo += str;
	}
	if(bfullalign){
		while(is1>0 || is2>0){
			char c1 = '-', c2 = '-';
			if(is1 > 0) c1 = rnam1_std[pa->resid[--is1]];
			if(is2 > 0) c2 = rnam1_std[pb->resid[--is2]];
			sa1 = c1 + sa1;
			sa2 = c2 + sa2;
			sdis = ' ' + sdis;
		}
		while(ie1<nA-1 || ie2<nB-1){
			char c1 = '-', c2 = '-';
			if(ie1 < nA-1) c1 = rnam1_std[pa->resid[++ie1]];
			if(is2 < nB-1) c2 = rnam1_std[pb->resid[++ie2]];
			sa1 = sa1 + c1;
			sa2 = sa2 + c2;
			sdis = sdis + ' ';
		}
		sprintf(str, "%s\n", sa1.c_str()); sinfo += str;
		sprintf(str, "%s\n", sdis.c_str()); sinfo += str;
		sprintf(str, "%s\n", sa2.c_str()); sinfo += str;
	} else {
		sprintf(str, "%-4d %s %d\n", is1, sa1.c_str(), ie1); sinfo += str;
		sprintf(str, "%-4s %s\n", "", sdis.c_str()); sinfo += str;
		sprintf(str, "%-4d %s %d\n", is2, sa2.c_str(), ie2); sinfo += str;
	}
	}
}
//
int Salign::fit1(vector<int> *iali2){
	double xc1[3], xc2[3];
	int nfit = 0;
	for(int i=0; i<iali2[0].size(); i++){
		int i1=iali2[0][i], i2=iali2[1][i];
		if(i1<0 || i2<0) continue;
		wfit[nfit] = 1.;
		for(int m=0; m<3; m++){
			xn1[nfit][m] = (*xap)[i1][m];
			xn2[nfit][m] = (*xbp)[i2][m];
		}
		nfit ++;
	}
	if(nfit < 3) die("too small alignment size: %d", nfit);
	double ds = 0.;
	for(int i=0; i<nfit; i++) ds += wfit[i];
	if(ds > 1.0e-8){
		for(int i=0; i<nfit; i++) wfit[i] /= ds;
	} else {
		for(int i=0; i<nfit; i++) wfit[i] = 1./nfit;
	}
	translate2(nfit, wfit, xn1, xn1, xc1);
	translate2(nfit, wfit, xn2, xn2, xc2);
	quatfit(nfit, wfit, xn2, xn1, u);
	for(int m=0; m<3; m++) u[3][m] = xc2[m] - dot_product(u[m], xc1);
	rotmol(nA, u, (*xap), xn1);
	return 0;
}
void Salign::scoreOnly(){
	Rscore_max = Rscore_sv = -100;
	for(int m=0; m<2; m++) ialign[m].clear();
	assert(nA <= nB);
	if(bscoreOnly == 1){
		map<string, int> resnum2id;
		for(int i=0; i<pa->resnum.size(); i++){
			resnum2id[pa->resnum[i]] = i;
		}
		for(int i=0; i<pb->resnum.size(); i++){
			if(resnum2id.count(pb->resnum[i]) <= 0) continue;
			ialign[0].push_back(resnum2id[pb->resnum[i]]);
			ialign[1].push_back(i);
		}
	} else {
		for(int i=0; i<nA; i++){
			ialign[0].push_back(i);
			ialign[1].push_back(i);
		}
	}
	int nfit0=ialign[0].size(),  nfit=nfit0;
	if(nfit0 < 3) die("too small number of match residues for score: %d", nfit0);
	fit1(ialign);
	double rms=0; calrms(rms);
	if(iprint>=99) printf("initial RMSD: %.2f\n", rms);
//
	vector<int> iali2[2];
	for(int k=0; k<10; k++){
		for(int i=0; i<nfit0-nfit+1; i++){
			for(int m=0; m<2; m++) iali2[m].clear();
			for(int m=0; m<nfit; m++){
				iali2[0].push_back( ialign[0][i+m] );
				iali2[1].push_back( ialign[1][i+m] );
			}
			fit1(iali2);
			optimize_score(1);
		}
		if(nfit <= 4) break;
		nfit /= 4;
		if(nfit < 4) nfit = 4;
	}
	Rscore = Rscore_max;
}
