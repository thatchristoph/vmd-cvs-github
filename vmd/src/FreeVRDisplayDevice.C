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
 *      $RCSfile: FreeVRDisplayDevice.C,v $
 *      $Author: johns $        $Locker:  $                $State: Exp $
 *      $Revision: 1.30 $      $Date: 2011/06/10 05:04:24 $
 *
 ***************************************************************************
 * DESCRIPTION:
 *
 * a FreeVR specific display device for VMD
 ***************************************************************************/

#include <freevr.h> // include FreeVR library prototypes
#include "Inform.h"
#include "FreeVRDisplayDevice.h"

// static string storage used for returning stereo modes
static const char *freevrStereoNameStr[1] = {"FreeVR"};

///////////////////////////////  constructor
FreeVRDisplayDevice::FreeVRDisplayDevice(void) : OpenGLRenderer("FreeVR") {
  stereoNames = freevrStereoNameStr;
  stereoModes = 1;
  doneGLInit = FALSE;    
  num_display_processes  = vrContext->config->num_windows;

  // XXX migrated some initialization code into the constructor due
  //     to the order of initialization in CAVE/FreeVR builds
  aaAvailable = TRUE;               // enable antialiasing
  cueingAvailable = FALSE;          // disable depth cueing
  cullingAvailable = FALSE;         // disable culling
  ext->hasstereo = TRUE;            // stereo is on initially
  ext->stereodrawforced = FALSE;    // no need for force stereo draws

  ogl_useblendedtrans = 0;
  ogl_transpass = 0;
  ogl_useglslshader = 0;
  ogl_acrobat3dcapture = 0;
  ogl_lightingenabled = 0;
  ogl_rendstateserial = 1;    // force GLSL update on 1st pass
  ogl_glslserial = 0;         // force GLSL update on 1st pass
  ogl_glsltoggle = 1;         // force GLSL update on 1st pass
  ogl_glslmaterialindex = -1; // force GLSL update on 1st pass
  ogl_glslprojectionmode = DisplayDevice::PERSPECTIVE;
  ogl_glsltexturemode = 0;    // initialize GLSL projection to off


  // leave everything else up to the freevr_gl_init_fn
}

///////////////////////////////  destructor
FreeVRDisplayDevice::~FreeVRDisplayDevice(void) {
  // nothing to do
}


/////////////////////////////  public routines  //////////////////////////

// set up the graphics on the seperate FreeVR displays
void FreeVRDisplayDevice::freevr_gl_init_fn(void) {
  setup_initial_opengl_state();     // do all OpenGL setup/initialization now

  // follow up with mode settings
  aaAvailable = TRUE;               // enable antialiasing
  cueingAvailable = FALSE;          // disable depth cueing
  cullingAvailable = FALSE;         // disable culling
  ext->hasstereo = TRUE;            // stereo is on initially
  ext->stereodrawforced = FALSE;    // no need for force stereo draws

  glClearColor(0.0, 0.0, 0.0, 0.0); // set clear color to black

  aa_on();                          // force antialiasing on if possible
  cueing_off();                     // force depth cueing off

  // set default settings
  set_sphere_mode(sphereMode);
  set_sphere_res(sphereRes);
  set_line_width(lineWidth);
  set_line_style(lineStyle);

  clear();                          // clear screen
  update();                         // swap buffers

  // we want the CAVE to be centered at the origin, and in the range -1, +1
  (transMat.top()).translate(0.0, 3.0, -2.0);
  (transMat.top()).scale(VMD_PI);

  doneGLInit = TRUE;                // only do this once
}

void FreeVRDisplayDevice::set_stereo_mode(int) {
  // cannot change to stereo mode in FreeVR, it is setup at init time
}

void FreeVRDisplayDevice::normal(void) {
  // prevent the OpenGLRenderer implementation of this routine
  // from overriding the projection matrices provided by the
  // FreeVR library.
}

// special render routine to check for graphics initialization
void FreeVRDisplayDevice::render(const VMDDisplayList *cmdlist) {
  if(!doneGLInit) {
    freevr_gl_init_fn();
  }

  // prepare for rendering
  glPushMatrix();
  multmatrix((transMat.top()));  // add our FreeVR adjustment transformation

  // update the cached transformation matrices for use in text display, etc.
  // In FreeVR, we have to do this separately for all of the processors.
  // Would be nice to do this outside of the render routine however,
  // amortized over several Displayables.
  glGetFloatv(GL_PROJECTION_MATRIX, ogl_pmatrix);
  glGetFloatv(GL_MODELVIEW_MATRIX, ogl_mvmatrix);
  ogl_textMat.identity();
  ogl_textMat.multmatrix(ogl_pmatrix);
  ogl_textMat.multmatrix(ogl_mvmatrix);

  // call OpenGLRenderer to do the rest of the rendering the normal way
  OpenGLRenderer::render(cmdlist);
  glPopMatrix();
}

// update after drawing
void FreeVRDisplayDevice::update(int do_update) {
  // XXX don't do buffer swaps in FreeVR!!!
  //     Though not well documented, it is implicitly illegal 
  //     to call glxSwapBuffers() or to call glDrawBuffer() 
  //     in a FreeVR application, since FreeVR does this for you.
}

