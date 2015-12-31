#include "protein2.h"

namespace Salign_PARAMS{
// alpha is the normalized factor(as seen in paper); D00 is "D0"; Denv the distance cutoff used to determine the "environment residues"
	double Alpha=0.3, D00=4., Denv=12., cutoff=-1;
// score_type decids the type of calculated alignment scores (SP, TM, GDT-scores are supported)
// fragsize is the length of fragment for the original alignment trials during finding the best alignment
	int score_type=iSP, fragsize=20; 
// iprint controls the details to print, the bigger the more printed
	int iprint=-9999;
// bscoreOnly: 0, structure alignment; others, score only (using predefined alignment)
// 1,scored according to resi No.; 2,scored according to residues sequentially
	int bscoreOnly=0;
	bool rmsOnly=0;
	bool bfullalign=0, nonSeqAlign=0;
// repAtom defines the representative atom from input structure to use for alignment
    string repAtom = "CA";
}
namespace PARAMS{
	using namespace Salign_PARAMS;
	vector<string> Tlist, Qlist;
	string fali, idir, odir;
	int bpairlist=0, bcheck=0;
	vector<string> folds;
}
static string runtype = "";
using namespace PARAMS;
int DEBUG = 0;
void rdlist(string sdir, string flist, vector<string> &slist){
	string line; vector<string> ss;
	ifstream fs(flist.c_str());
	if(sdir.size()>0 && sdir[sdir.size()-1]!='/') sdir += '/';
	while(getline(fs, line)){
		if(line[0] == '#') continue;
		int n = str2dat(line, ss);
		if(n < 1) continue;
		slist.push_back(sdir + ss[0]);
	}
	fs.close();
}
inline string getdir(string sdir0){
	string sdir = sdir0;
	if(sdir.size()>0 && sdir[sdir.size()-1]!='/') sdir += '/';
	return sdir;
}
void rdparams(int argc, char *argv[]){
	if(argc < 3) 
		die("Usage: RUN [-pair pdb1 pdb2|-pairlist sdir list] [-iprint 1] [-SP|-TM|-GDT] [-NS] [-fast] [-cutoff cut] [-repAtom atom]");

	Tlist.clear(); Qlist.clear();
	fali = idir = odir = "";
	int i0 = 1;
	while(i0 < argc){
		if(i0==1 && argv[i0][0] != '-'){
			Qlist.push_back(argv[i0 ++]);
			Tlist.push_back(argv[i0 ++]);
			iprint = 100;
			continue;
		}
		string opt1 = argv[i0];
		if(opt1 == "-tlist"){
			rdlist(argv[i0+1], argv[i0+2], Tlist);
			i0 += 2;
		}else if(opt1 == "-qlist"){
			rdlist(argv[i0+1], argv[i0+2], Qlist);
			i0 += 2;
		}else if(opt1 == "-pairlist"){
			bpairlist = 1;
			rdlist(argv[i0+1], argv[i0+2], Tlist);
			i0 += 2;
		}else if(opt1 == "-pair"){
			if(iprint == -9999) iprint = 100;
			Qlist.push_back(idir + argv[++ i0]);
			Tlist.push_back(idir + argv[++ i0]);
		}else if(opt1 == "-t"){
			for(; i0<argc-1; i0++){
				if(argv[i0+1][0] == '-') break;
				Tlist.push_back(argv[i0+1]);
			}
		} else if(opt1 == "-q"){
			for(; i0<argc-1; i0++){
				if(argv[i0+1][0] == '-') break;
				Qlist.push_back(argv[i0+1]);
			}
		}
		else if(opt1 == "-idir") idir = getdir(argv[++i0]);
		else if(opt1 == "-odir") odir = getdir(argv[++i0]);
		else if (opt1 == "-iprint") iprint = atoi(argv[++i0]);
		else if(opt1 == "-scoreOnly") bscoreOnly = 1;
		else if(opt1 == "-scoreOnly2") bscoreOnly = 2;
		else if(opt1 == "-rmsOnly") rmsOnly = true;
		else if (opt1 == "-a") Alpha = strtod(argv[++i0], NULL);
		else if (opt1 == "-d0") D00 = strtod(argv[++i0], NULL);
		else if(opt1 == "-frag") fragsize = atoi(argv[++i0]);
		else if(opt1 == "-fast") fragsize = 0;
		else if (opt1 == "-denv") Denv = strtod(argv[++i0], NULL);
//		else if (opt1 == "-outpdb") outpdb = argv[++i0];
		else if(opt1 == "-check") bcheck = atoi(argv[++i0]);
		else if(opt1 == "-fali") fali = argv[++i0];
		else if(opt1 == "-cutoff") cutoff = strtod(argv[++i0], NULL);
		else if(opt1 == "-plist") {fali = argv[++i0]; runtype="plist";}
		else if(opt1 == "-SP") score_type = iSP;
		else if(opt1 == "-GDT") score_type = iGDT;
		else if(opt1 == "-TM") score_type = iTM;
		else if(opt1 == "-LG") score_type = iLG;
		else if(opt1 == "-fullalign") bfullalign = 1;
		else if(opt1 == "-NS") nonSeqAlign = 1;
        else if(opt1 == "-repAtom") repAtom = argv[++i0];
		else die("unknown option: %s", argv[i0]);
		i0 ++;
	}
	if(iprint == -9999) iprint = 0;
	if(bpairlist) Qlist = Tlist;
//	fprintf(stderr, "protein number: %d %d\n", Qlist.size(), Tlist.size());
}
void align_list_mp(){
// read the list of "template" proteins
	vector<Protein2*> Tpro;
	for(int m=0; m<Tlist.size(); m++){
		Protein2 *p1 = new Protein2(Tlist[m]);
		Tpro.push_back(p1);
	}
// scan the list of "query" proteins
	FILE *fp = stdout;
	for(int j=Qlist.size()-1; j>=0; j--){
		string bn = basename((char*)Qlist[j].c_str());
		if(odir != ""){
			string fn = odir + bn + ".sp";
			if(file_existed(fn)) continue;
			fp = openfile(fn, "w");
		}
//
		Protein2 *p1;
		int ne = Tpro.size();
		if(bpairlist) {ne = j+1; p1 = Tpro[j];}
		else p1 = new Protein2(Qlist[j]);
		vector<string> sdim (ne, "");
#if defined MP
	#pragma omp parallel for schedule(dynamic, 1)
#endif
		for(int i=0; i<ne; i++){
//
// here is the real start to do the alignment
			Salign *sa1 = new Salign();
			sa1 -> init(p1, Tpro[i], 1);
			sa1 -> Run_all();
#if defined MP
			sa1 -> print_max(sdim[i]);
#else
			sa1 -> print_max(fp);
#endif
			delete sa1;
		}
#if defined MP
		for(int i=0; i<ne; i++){
			fprintf(fp, "%s", sdim[i].c_str());
		}
#endif
		if(odir!="" && fp!=stdout) {
			fprintf(fp, "END\n");
			fclose(fp);
		}
		if(! bpairlist) delete p1;
	}
}
void align_pair1(){
	vector<Protein2*> Tpro;
	Tpro.resize(Tlist.size(), NULL);
	FILE *fpo = stdout;
	Salign *sa1 = new Salign();
	for(int j=0; j<Qlist.size(); j++){
		Protein2 *p1 = new Protein2(Qlist[j]);
		int ne = Tpro.size();
		if(bpairlist) ne = j+1;
		for(int i=0; i<ne; i++){
			if(Tpro[i] == NULL) Tpro[i] = new Protein2(Tlist[i]);
			sa1 -> init(p1, Tpro[i], 1);
			sa1 -> Run_all();
			sa1 -> print_max(fpo);
		}
		delete p1;
	}
	delete sa1;
}
// only run for file containing pairs
void align_plist(string idir, string fali){
	map<string, Protein2*> Plist;
	FILE *fp = openfile(fali, "r");
	char str[501], p1[501], p2[501];
	Salign *sa1 = new Salign();
	while(fgets(str, 500, fp) != NULL){
		sscanf(str, "%s%s", p1, p2);
		string fn1 = idir + p1, fn2 = idir + p2;
		if(Plist.count(fn1) < 1) Plist[fn1] = new Protein2(fn1);
		if(Plist.count(fn2) < 1) Plist[fn2] = new Protein2(fn2);
		sa1 -> init(Plist[fn1], Plist[fn2], 1);
		sa1 -> Run_all();
		sa1 -> print_max(stdout);
	}
}
int main(int argc, char *argv[]){
	rdparams(argc, argv);
	if(Tlist.size() < 1) die("no tpl selected");
	if(Qlist.size() < 1) die("no query selected");
	if(runtype == "") {
// normal mode comparison
		align_list_mp();

	} else if(runtype == "plist"){
// the mode to run one list of pairs
// idir is the directory to put all proteins and fali includes two columns of pdb
//
		align_plist(idir, fali);
	} else die("unknown running mode: %s", runtype.c_str());
}
