// $Id: Isotropic3DElemDataCard.C,v 1.6 2006-09-05 20:41:51 manav Exp $

// FESystem includes
#include "Properties/Isotropic3DElemDataCard.h"


Isotropic3D_ElemDataCard::Isotropic3D_ElemDataCard():
IsotropicElemDataCard()
{
  
}



Isotropic3D_ElemDataCard::~Isotropic3D_ElemDataCard()
{
  
}




std::istream& operator>> (std::istream& input, Isotropic3D_ElemDataCard& solid_elem)
{
  return solid_elem.readFromInputStream(input);
}
  
