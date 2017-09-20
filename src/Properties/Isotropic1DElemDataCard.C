// $Id: Isotropic1DElemDataCard.C,v 1.6.6.1 2007-03-14 22:05:03 manav Exp $

// FESystem includes
#include "Properties/Isotropic1DElemDataCard.h"



Isotropic1D_ElemDataCard::Isotropic1D_ElemDataCard():
IsotropicElemDataCard()
{
  
}




Isotropic1D_ElemDataCard::~Isotropic1D_ElemDataCard()
{
  
}







std::istream& operator>> (std::istream& input, Isotropic1D_ElemDataCard& line_elem)
{
  return line_elem.readFromInputStream(input);
}


void
Isotropic1D_ElemDataCard::getFactor(DenseMatrix<double>& factor,
                                    const unsigned int factor_enum_ID)
{
  // params not used
  (void) factor;
  (void) factor_enum_ID;

  Assert(false, ExcPureFunctionCalled());
}



void
Isotropic1D_ElemDataCard::getFactorSensitivityForGlobalParameter
(DenseMatrix<double>& factor,
 const unsigned int factor_enum_ID,
 const unsigned int global_param_ID)
{
  // params not used
  (void) factor;
  (void) factor_enum_ID;
  (void) global_param_ID;

  Assert(false, ExcPureFunctionCalled());
}


void
Isotropic1D_ElemDataCard::getFactorSensitivityForLocalParameter
(DenseMatrix<double>& factor,
 const unsigned int factor_enum_ID,
 const unsigned int param_enum_ID)
{
  // params not used
  (void) factor;
  (void) factor_enum_ID;
  (void) param_enum_ID;

  Assert(false, ExcPureFunctionCalled());
}

