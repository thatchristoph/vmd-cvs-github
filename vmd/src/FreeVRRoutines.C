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
 *	$RCSfile: FreeVRRoutines.C,v $
 *	$Author: johns $	$Locker:  $		$State: Exp $
 *	$Revision: 1.21 $	$Date: 2011/06/11 03:24:16 $
 *
 ***************************************************************************
 * DESCRIPTION:
 *
 * routines to get memory from and return memory to the 
 * FreeVR shared memory arena
 ***************************************************************************/

#include <stdlib.h>
#include "Inform.h"
#include "FreeVRRoutines.h"
#include "VMDApp.h"
#include "FreeVRDisplayDevice.h"
#include "FreeVRScene.h"

#include <freevr.h>

void *malloc_from_FreeVR_memory(size_t size) {
  return vrShmemAlloc(size);
}

void free_to_FreeVR_memory(void *data) {
  vrShmemFree(data);
}

// get megs o' memory from FreeVR, and create the arena
// Warning:  Don't make me do this twice.
void grab_FreeVR_memory(size_t megs) {
  size_t size = ((megs>1) ? megs : 1) * 1024L * 1024L;
  size_t sz=0;

  while (((sz = vrShmemInit(size)) == 0) && (size > 64*1024*1024)) {
    msgErr << "Failed to create FreeVR arena of size " 
           << (size / (1024*1024)) 
           << ", reducing allocation by half." << sendmsg;
    size >>= 1; // cut allocation in half
  }
 
  if (sz == 0) 
    msgErr << "Failed to create FreeVR arena.  We're gonna die!" << sendmsg;
  else
    msgInfo << "Created arena, size " << (sz / (1024*1024)) 
            << "MB." << sendmsg;
}


// set up the graphics, called from FreeVRInitApplication
void freevr_gl_init_fn(void) {
}

static Scene *freevrscene;
static DisplayDevice *freevrdisplay;

void set_freevr_pointers(Scene *scene, DisplayDevice *display) {
  freevrscene = scene;
  freevrdisplay = display;
}

// call the child display renderer, and wait until they are done
void freevr_renderer(void) {
  freevrscene->draw(freevrdisplay);
}

