// $Id: DesignParameter.C,v 1.3 2006-09-05 20:41:41 manav Exp $

// FESystem includes
#include "DesignData/DesignParameter.h"


DesignData::DesignParameter::DesignParameter():
ID(FESystemNumbers::InvalidID),
name(),
value(0.0),
sensitivity_method(FESystemNumbers::InvalidID),
finite_difference_step_size(0.0)
{}


DesignData::DesignParameter::~DesignParameter()
{}


