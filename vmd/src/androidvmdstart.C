/***************************************************************************
 *cr
 *cr            (C) Copyright 1995-2011 The Board of Trustees of the
 *cr                        University of Illinois
 *cr                         All Rights Reserved
 *cr
 ***************************************************************************/

/***************************************************************************
 * RCS INFORMATION:
 *
 *      $RCSfile: androidvmdstart.C,v $
 *      $Author: johns $        $Locker:  $             $State: Exp $
 *      $Revision: 1.3 $      $Date: 2012/10/18 15:23:11 $
 *
 ***************************************************************************
 * DESCRIPTION:
 *  Android startup code
 ***************************************************************************/

// only compile this file if we're building on Android
#if defined(ANDROID)
#include "androidvmdstart.h"
#include "vmd.h"

//
// VMD is wrapped in a JNI shared object and called by Java...
//
#include <stdio.h>
#include <jni.h>

#if !defined(VMD_JNI_CLASSNAME)
#define VMD_JNI_CLASSNAME "edu/uiuc/VMD/VMD"
#endif
#if !defined(VMD_JNI_WRAPFUNC)
#define VMD_JNI_WRAPFUNC Java_edu_uiuc_VMD_VMD_nativeMain
#endif


/* ------------------------------------------------------------ */
/* output is allocated elsewhere */
/* key: FilesDir    - gets internal phone memory space for the app */
/* key: VMDDIR      - gets internal phone memory space, where VMD files live*/
/* key: ExternalStorageDir  - gets SDcard memory space (global, use with care)*/
char *getPlatformValue(JNIEnv *env, jobject thiz, const char *key, char *output) {

   jclass clazz = (env)->FindClass("edu/uiuc/VMD/VMD");

   /* find the method that takes a java String and returns a String */
   /* this call is expensive */
   jmethodID getPlatformValue = (env)->GetMethodID(clazz,
                                     "getPlatformValue",
                                     "(Ljava/lang/String;)Ljava/lang/String;");

   /* NewStringUTF makes a java string from a char* */
   jobject j = (env)->CallObjectMethod(thiz, getPlatformValue,
                                       (env)->NewStringUTF(key));
   const char* str = (env)->GetStringUTFChars((jstring) j, NULL);

   strcpy(output, str);

   (env)->ReleaseStringUTFChars((jstring) j, str);
   return output;
}



void logtojava(JNIEnv *env, jobject thiz, const char *logstring) {
  // c++ version JNI interface:
  jclass clazz = (env)->FindClass(VMD_JNI_CLASSNAME);
  jmethodID logOutput = (env)->GetMethodID(clazz, "logOutput",
                                           "(Ljava/lang/String;)V");

  /* this actually logs a char*, logstring */
  (env)->CallVoidMethod(thiz, logOutput, (env)->NewStringUTF(logstring));
}



// XXX gross disgusting hack!!!!
JNIEnv *global_jnienv; // XXX hack! 
jobject global_thiz;   // XXX hack! 

extern "C" {

//
// Wrapper function to hide use of the cached global state until
// until we make appropriate changes so that the JNI launcher has
// a mechanism to provide VMD with the JNI objects for use by calls
// back to android APIs.
//
void log_android(const char *prompt, const char * msg) {
  char logstring[2048];

  strncpy(logstring, prompt, sizeof(logstring)-2);
  strcat(logstring, msg);
  strcat(logstring, "\n");

  logtojava(global_jnienv, global_thiz, logstring);
}


//
// This is the main JNI wrapper function.
// Contains startup code, neverending loop, shutdown code, etc...
//
void VMD_JNI_WRAPFUNC(JNIEnv* env, jobject thiz) {
  char* rargv[10];
 
  global_jnienv = env; // XXX this is a hack!
  global_thiz = thiz;  // XXX this is a hack!

  fprintf(stderr, "--stderr fprintf---------------------------------\n");
  printf("---regular printf----------------------------\n");
  fflush(stdout);
  log_android("", "--Log event ---------------------");

#if 1
  printf("VMD Android platform info:\n");
  printf("  sizeof(char): %d\n", sizeof(char));
  printf("  sizeof(int): %d\n", sizeof(int));
  printf("  sizeof(long): %d\n", sizeof(long));
  printf("  sizeof(void*): %d\n", sizeof(void*));
  fflush(stdout);
#endif

  char tmp[8192];
  const char * vmddir = NULL;

  // set to a worst-case guess until we have something better.
  vmddir = "/data/data/edu.uiuc.VMD/files/vmd";

#if 1
  // Query Android for app directories and files here...
  char androidappdatadir[8192];

  memset(androidappdatadir, 0, sizeof(androidappdatadir));
  getPlatformValue(global_jnienv, global_thiz, "FilesDir", androidappdatadir);

  if (strlen(androidappdatadir) > 0) {
//    log_android("ANDROID APP DIR: ", androidappdatadir);
    strcat(androidappdatadir, "/vmd");
    vmddir = androidappdatadir;
  }
#endif

  if (vmddir == NULL) {
    return; // fail/exit
  }

  if (!getenv("VMDDIR")) {
    setenv("VMDDIR", vmddir, 1);
  }

  if (!getenv("TCL_LIBRARY")) {
    strcpy(tmp, vmddir);
    strcat(tmp, "/scripts/tcl");
    setenv("TCL_LIBRARY", tmp, 1);
  }

  if (!getenv("TK_LIBRARY")) {
    strcpy(tmp, vmddir);
    strcat(tmp, "/scripts/tk");
    setenv("TK_LIBRARY", tmp, 1);
  }

  if (!getenv("PYTHONPATH")) {
    strcpy(tmp, vmddir);
    strcat(tmp, "/scripts/python");
    setenv("PYTHONPATH", tmp, 1);
  } else {
    strcpy(tmp, getenv("PYTHONPATH"));
    strcat(tmp, ":");
    strcat(tmp, vmddir);
    strcat(tmp, "/scripts/python");
    setenv("PYTHONPATH", tmp, 1);
  }

  if (!getenv("STRIDE_BIN")) {
    strcpy(tmp, vmddir);
#if defined(ARCH_ANDROIDARMV7A)
    strcat(tmp, "/stride_ANDROIDARMV7A");
#else
#error unhandled compilation scenario
#endif
    setenv("STRIDE_BIN", tmp, 1);
  }

  if (!getenv("SURF_BIN")) {
    strcpy(tmp, vmddir);
#if defined(ARCH_ANDROIDARMV7A)
    strcat(tmp, "/surf_ANDROIDARMV7A");
#else
#error unhandled compilation scenario
#endif
    setenv("SURF_BIN", tmp, 1);
  }

  if (!getenv("TACHYON_BIN")) {
    strcpy(tmp, vmddir);
#if defined(ARCH_ANDROIDARMV7A)
    strcat(tmp, "/tachyon_ANDROIDARMV7A");
#else
#error unhandled compilation scenario
#endif
    setenv("TACHYON_BIN", tmp, 1);
  }

  rargv[0] = "VMD.so";
#if 1
  rargv[1] = "1e79";
#elif 1
  rargv[1] = "/data/data/edu.uiuc.VMD/files/alanin.pdb";
#else
  rargv[1] = "-h";
#endif
  rargv[2] = NULL;

  VMDmain(2, rargv); /* launch VMD... */

  log_android("", "--Log event ---------------------");
  fprintf(stderr, "--stderr fprintf---------------------------------\n");
  printf("---regular printf----------------------------\n");
  fflush(stdout);
}

} // extern "C"

#endif

