#include "sim2spImpl.h"

double SPalignNS_RMSD(JNIEnv *env, jobject thisObj, jstring pdb1, jstring pdb2){

	const char *s1 = env->GetStringUTFChars(pdb1, 0);
	const char *s2 = env->GetStringUTFChars(pdb2, 0);
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
