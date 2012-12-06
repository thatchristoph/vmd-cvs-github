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
 *	$RCSfile: P_FreeVRTracker.h,v $
 *	$Author: johns $	$Locker:  $		$State: Exp $
 *	$Revision: 1.15 $	$Date: 2010/12/16 04:08:30 $
 *
 ***************************************************************************
 * DESCRIPTION:
 * This is Paul's new Tracker code -- pgrayson@ks.uiuc.edu
 * Representation of the Tracker in FreeVR
 ***************************************************************************/

/// VMDTracker subclass that interfaces to the FreeVR wand
class FreeVRTracker : public VMDTracker {
 public:
  inline FreeVRTracker() : VMDTracker() {}; // FreeVR needs no initialization
  virtual VMDTracker *clone() { return new FreeVRTracker; }
  const char *device_name() const { return "freevrtracker"; }
  virtual void update();
  inline virtual int alive() { return 1; }

 protected:
  virtual int do_start(const SensorConfig *);
};

