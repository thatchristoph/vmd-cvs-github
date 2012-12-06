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
 *      $RCSfile: vmd.h,v $
 *      $Author: johns $        $Locker:  $             $State: Exp $
 *      $Revision: 1.9 $       $Date: 2012/10/17 04:07:15 $
 *
 ***************************************************************************
 * DESCRIPTION:
 *
 * Main program interface header.
 *
 ***************************************************************************/

#ifndef VMD_H_
#define VMD_H_

// This header file contains all the functions exported by the VMD library
// for various VMD implementations; e.g. standalone executable, Python 
// module, etc.

#ifdef ANDROID
// prototype for main program JNI shared library entry point 
int VMDmain(int argc, char **argv);
#endif

#include "VMDApp.h"

// Perform various one-time startup operations prior to creating a
// VMDApp instance.  Reads the command line arguments, processes environnment
// variables, and prepares a set startup operations that will be carried
// out by VMDreadStartup.  Upon return, *argc will be set to the number of 
// unprocessed arguments that should be passed to VMDinit, including argv[0]
// (the path to the executable) and argv will hold the arguments themselves.
extern int VMDinitialize(int *argc, char ***argv);

// Return the display type name that was chosen by VMDinitialize()
extern const char *VMDgetDisplayTypeName();

// Return the display location and display size that were chosen by
// VMDinitialize()
extern void VMDgetDisplayFrame(int *loc, int *size);

// read various default settings defined in scripts/vmd/.  Should be done
// after a successful VMDinitialize, but before VMDreadStartup().
void VMDreadInit(VMDApp *app);
 
// Get the unprocessed command line arguments 
// Read startup files (.vmdrc, -e files), and load molecules specified on
// the command line.  Must be called _after_ calling VMDinitialize()
extern void VMDreadStartup(VMDApp *);

// Shut down VMD properly, after all VMDApp instances have been deleted
extern void VMDshutdown();

// check for Fltk events, if available.
extern void VMDupdateFltk();

#endif
