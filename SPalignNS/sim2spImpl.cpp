#include "sim2spImpl.h"

double SPalignNS_RMSD(JNIEnv *env, jobject thisObj, jstring pdb1, jstring pdb2){

    const char *s1 = env->GetStringUTFChars(pdb1, 0);
    const char *s2 = env->GetStringUTFChars(pdb2, 0);

    //std::istringstream iss(s1);
    //std::string line;    
		//double res = 0;
    //while (std::getline(iss,line)) {
    //    //std::cout << line << std::endl;
		//		res++;
    //}
		//return res;
		
    std::string result;    
		Protein2 *p1 = new Protein2(s1);
		Protein2 *p2 = new Protein2(s2);
		Salign *sa1 = new Salign();
		sa1 -> init(p1, p2, 1);
		sa1 -> Run_all();
		sa1 -> print_max(result);
		env->ReleaseStringUTFChars(pdb1, s1);
		env->ReleaseStringUTFChars(pdb2, s2);
		return ::atof(result.c_str());
}
//int main(int argc, char *argv[]){
//	std::cout << "in\n";
//  std::string content;    
//  std::string result;    
//	std::ifstream ifs(argv[0]);
//	ifs.seekg(0, std::ios::end);   
//  content.reserve(ifs.tellg());
//  ifs.seekg(0, std::ios::beg);
//  content.assign((std::istreambuf_iterator<char>(ifs)),
//                std::istreambuf_iterator<char>());
//	Protein2 *p1 = new Protein2(content);
//	std::cout << "end p1\n";
//	std::ifstream ifs2(argv[0]);
//	ifs2.seekg(0, std::ios::end);   
//  content.reserve(ifs2.tellg());
//  ifs2.seekg(0, std::ios::beg);
//  content.assign((std::istreambuf_iterator<char>(ifs2)),
//                std::istreambuf_iterator<char>());
//	Protein2 *p2 = new Protein2(content);
//	std::cout << "end p2\n";
//	Salign *sa1 = new Salign();
//	sa1 -> init(p1, p2, 1);
//	sa1 -> Run_all();
//	sa1 -> print_max(result);
//	std::cout << "e " << ::atof(result.c_str())<< "\n";
//}
