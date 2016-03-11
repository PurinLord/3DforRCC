#include <jni.h>

#include "rccto3d_optimisation_SimulatedAnnealing3DProtFromRCC.h"
#include "sim2spImpl.h"

JNIEXPORT jdouble JNICALL Java_rccto3d_optimisation_SimulatedAnnealing3DProtFromRCC_calcRMSD_1New_1fast
  (JNIEnv *env, jobject thisObj, jstring pdb1, jstring pdb2){
		return SPalignNS_RMSD(env, thisObj, pdb1, pdb2);
	}
