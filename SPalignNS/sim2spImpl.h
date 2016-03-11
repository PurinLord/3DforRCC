#ifndef __sim2spImpl__
#define __sim2spImpl__ 

#include <jni.h>

#include <stdio.h>
#include <string>
#include <sstream>
#include <iostream>

#include "protein2.h"

extern "C"{
	double SPalignNS_RMSD(JNIEnv *env, jobject thisObj, jstring pdb1, jstring pdb2);
}
#endif
