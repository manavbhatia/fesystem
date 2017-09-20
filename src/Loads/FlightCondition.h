// $Id:$
/*
 *  FlightCondition.h
 *  fesystem_xcode
 *
 *  Created by Manav Bhatia on 12/25/08.
 *  Copyright 2008 Virginia Tech. All rights reserved.
 *
 */

#ifndef __fesystem_flight_condition_h__
#define __fesystem_flight_condition_h__

// libmesh includes
#include "geom/point.h"


namespace Loads {
  class FlightCondition {

  public:
    
    /// contructor
    FlightCondition();

    /// copy contructor
    FlightCondition(const Loads::FlightCondition& flt_cond);

    /// destructor 
    ~FlightCondition();
    
    /// @returns the flight dynamic pressure
    double getDynamicPressure() const;
    
    /// @returns the flight mach number
    double getMachNumber() const;

    /// @returns the flight density
    double getDensity() const;
    
    /// @returns the speed of sound at flight condition
    double getSpeedOfSound() const;

    /// @returns the fluid flow vector
    const Point& getFluidFlowVector() const;

    /// sets the value of the flight dynamic vector
    void setDynamicPressure(const double dyn_press);
    
    /// sets the flight mach number
    void setMachNumber(const double m_num);
    
    /// sets the flight density
    void setDensity(const double dens);

    /// sets the fluid flow vector
    void setFluidFlowVector(const Point& vec);
    
    
  protected:
    
    /// dynamic pressure
    double dynamic_pressure;
    
    /// mach number
    double mach_number;
    
    /// density
    double density;
    
    /// fluid flow vector
    Point fluid_flow_vector;
    
  };
}



#endif // __fesystem_flight_condition_h__


