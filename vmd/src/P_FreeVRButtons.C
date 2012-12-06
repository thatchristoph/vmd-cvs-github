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
 *	$RCSfile: P_FreeVRButtons.C,v $
 *	$Author: johns $	$Locker:  $		$State: Exp $
 *	$Revision: 1.15 $	$Date: 2011/06/08 21:18:10 $
 *
 ***************************************************************************
 * DESCRIPTION:
 * This is Paul's new Tracker code -- pgrayson@ks.uiuc.edu
 ***************************************************************************/

#include <freevr.h>
#include "P_FreeVRButtons.h"

int FreeVRButtons::do_start(const SensorConfig *) {
  // XXX Somehow check that a FreeVR environment exists.  If it doesn't,
  // return false.
  if ( 0 ) {
    return 0;  // no FreeVR, cannot run FreeVR buttons.
  }
  return 1;    // FreeVR is active.
}

void FreeVRButtons::update() {
  // the mapping of buttons is different in FreeVR than CAVElib
  stat[0]=vrGet2switchValue(1);
  stat[1]=vrGet2switchValue(2);
  stat[2]=vrGet2switchValue(3);
  stat[3]=vrGet2switchValue(0);
}

