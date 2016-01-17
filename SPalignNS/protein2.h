#ifndef _PROTEIN
#define _PROTEIN

#include "sp_misc.h"
#include "sp_type.h"

class Protein2{
	string name;
	int nres;
	vector<int> resid;
	vector<string> resnum;
	vector<Xvec> x;
	vector<int> ssec;
public:
	string info_err;
	Protein2(string fn){rdpdb(fn);};
	Protein2(int nr, string, double [][3]);
	void rdpdb(string fn);
	void wrpdb(string fn);
	void wrpdb(FILE*);
	int getnres(){return nres;}
	double distance2(int,int);
	void calSS();
	void calineib(double r2);
	void calnneib0(vector<int> &idx);
	void calnneib(vector<int> &idx, double);
	string getname(){return name;}
	~Protein2(){};
	friend class Salign;
	friend class NSalign;
};
void quatfit(int n, vector<double> &w, vector<Xvec> &x1, vector<Xvec> &x2, double u[][3]);
void translate(int na, vector<Xvec> &x0, vector<Xvec> &xn, double *xc);
void translate2(int na, vector<double> wfit, vector<Xvec> &x0, vector<Xvec> &xn, double *xc);
void rotmol(int na, double u[][3], vector<Xvec> &x0, vector<Xvec> &xn);
inline double distance2(double *xap, double *xbp);
//
//
class Salign{
	int ma, mb;
	Protein2 *pa, *pb;
	bool breverse;
	int nA, nB, Lmin;
	vector<Xvec> *xap, *xbp;		// quote for abbrev.
//
	double D0, D0_2;				// Rcut^2 from len_ali
	vector<Xvec> xn1, xn2;		// for changed coords
	vector<double> wfit;
// for DP
	double gap0, gap1;
	double Rscore, u[4][3];		//	current score; rotation matrix
	double scores_all[20];		// series of scores saved as ieLA... below
	double SO;
	vector< vector<short> > Idir;					// direction matrix for DP
	vector< vector<double> > Smat, Rmat;		// matrix for DP, sum & unit value
	vector<int> ialign[3];		// the aligned idx, only 2 are used 
// saved & max
	double Rscore_sv, u_sv[4][3];
	vector<int> ialign_sv[3];
	double Rscore_max, u_max[4][3];
	vector<int> ialign_max[3];
public:
	Salign():ma(0),mb(0){}
//	Salign(Protein2 *a, Protein2*b):ma(0),mb(0){init(a,b);}
	~Salign(){delete_x();}
	void init(Protein2 *a, Protein2*b, int brev);
	void init(Protein2 *a, Protein2*b){init(a,b,0);}
	void initalign(string smat[], string sali[]);
	void delete_x();
	void alloc_x();
	double optimize_score(int,int);
	double optimize_align();
	void Run_all();
	void Run1();
	void Run2();
	void Run3();
	void Run4();
	void Run_custom();
    void Run_NonSeq();
    void Refine_NonSeq();
	double run_Align1();
	double Align_DP();
	int calrms(double &rms, double cut2);
	int redo_fit(int flag, double c);
	int fit1(vector<int> *iali2);
	void calRmat();
	double calRscore();
	double getscore(){return Rscore;}
	void calscores_all();
	void calLali2(int*, double);
	bool save_align();
	void score_final();
	void print_max(FILE*);
	void print_max(string&);
	void prtali(FILE*);
	void prtali(string&);
	void checkAli1(string);
	void checkAlign(string idir, string fn, FILE *fpo, vector<Protein2*>&, Protein2 *p0);
	void pro_swap(int);
	void scoreOnly();
	void calSPscores(double*);

	friend class NSalign;
};
enum {iSP, iTM, iLG, iGDT};
enum {ieSP, ieTMa, ieLG, ieGDT, ieLA, ieRMS, ieLE, ieTMb, ieTMc, ieSEQ};
//
//                 TMc, TMa, TMb if iTM
//
//
static const int nres_std = 20;
static const char rnam3_std[][4] = 
		{"ALA","CYS","ASP","GLU","PHE","GLY","HIS","ILE","LYS","LEU",
		"MET","ASN","PRO","GLN","ARG","SER","THR","VAL","TRP","TYR"};
static const char rnam1_std0[] = "XACDEFGHIKLMNPQRSTVWY";
static const char *rnam1_std = rnam1_std0 + 1;
inline int aaDefine(string rn0, int DEBUG=0){
	string rn = rn0;
	if(rn=="HSD" || rn=="HSE") rn = "HIS";
	for(int i=0; i<nres_std; i++){
		if(rn==rnam3_std[i]) return i;
	}
	if(DEBUG > 0){
		fprintf(stderr, "Unrecognized residue name: %s\n", rn.c_str());
	}
	return -1;
}
inline int aaDefine1(char rn, int DEBUG=0){
	for(int i=0; i<nres_std; i++){
		if(rn==rnam1_std[i]) return i;
	}
	if(DEBUG > 0){
		fprintf(stderr, "Unrecognized residue name: %c\n", rn);
	}
	return -1;
}
//
namespace Salign_PARAMS{
	extern int iprint, inorm, score_type, fragsize, bscoreOnly;
	extern double Alpha, Beta, D00, Denv, cutoff;
	extern bool bfullalign, nonSeqAlign, rmsOnly, soOnly;
	extern string outpdb, repAtom;
}
using namespace Salign_PARAMS;

namespace PARAMS{
	extern vector<string> Tlist, Qlist;
	extern string fali, idir, odir;
	extern int bpairlist, bcheck;
	extern vector<string> folds;
}
using namespace PARAMS;

#endif
