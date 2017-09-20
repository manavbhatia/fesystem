// $Id:$
/*
 *  FlightCondition.cpp
 *  fesystem_xcode
 *
 *  Created by Manav Bhatia on 12/25/08.
 *  Copyright 2008 Virginia Tech. All rights reserved.
 *
 */

// FESystem includes
#include "FlightCondition.h"




Loads::FlightCondition::FlightCondition():
dynamic_pressure(0.0),
mach_number(0.0),
density(0.0)
{
  
}


Loads::FlightCondition::FlightCondition(const Loads::FlightCondition& flt_cond):
dynamic_pressure(flt_cond.dynamic_pressure),
mach_number(flt_cond.mach_number),
density(flt_cond.density),
fluid_flow_vector(flt_cond.fluid_flow_vector)
{
  
}


Loads::FlightCondition::~FlightCondition()
{
  
}



double
Loads::FlightCondition::getDynamicPressure() const
{
  return this->dynamic_pressure;
}



double
Loads::FlightCondition::getMachNumber() const
{
  return this->mach_number;
}


double
Loads::FlightCondition::getDensity() const
{
  return this->density;
}


double
Loads::FlightCondition::getSpeedOfSound() const
{
  return (sqrt(2.0 * this->dynamic_pressure / this->density) / this->mach_number); 
}




const Point& 
Loads::FlightCondition::getFluidFlowVector() const
{
  return this->fluid_flow_vector;
}



void
Loads::FlightCondition::setDynamicPressure(const double dyn_press)
{
  this->dynamic_pressure = dyn_press;
}



void
Loads::FlightCondition::setMachNumber(const double m_num)
{
  this->mach_number = m_num;
}


void
Loads::FlightCondition::setDensity(const double dens)
{
  this->density = dens;
}


void
Loads::FlightCondition::setFluidFlowVector(const Point& vec)
{
  this->fluid_flow_vector = vec;
}


