  /// \file Tracker.hh
/*
 *
 * Tracker.hh header template generated by fclass
 * Creation date : mer. juin 12 2013
 * Copyright (c) CNRS , IPNL
 *
 * All Right Reserved.
 * Use and copying of these libraries and preparation of derivative works
 * based upon these libraries are permitted. Any copy of these libraries
 * must include this copyright notice.
 * 
 * @author : rete
 */


#ifndef TRACKER_HH
#define TRACKER_HH

#include "Detector.hh"

#include <iostream> 
#include <string> 
#include <cstdlib> 
#include <cmath> 
#include <vector> 

namespace baboon {

/* 
 * Class Tracker
 * Inherits from base class Detector
 */ 

class Tracker : public Detector {

  public:
    /*! Default Constructor */
    Tracker();
    /*! Default Destructor */
    virtual ~Tracker();

  protected:


};  // class 

}  // namespace 

#endif  //  TRACKER_HH
