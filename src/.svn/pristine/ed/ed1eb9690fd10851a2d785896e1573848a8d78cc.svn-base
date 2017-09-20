// $Id: ElemPostProcessQty.C,v 1.7 2006-09-05 20:41:50 manav Exp $

// C++ includes
#include <string>
#include <sstream>
#include <cassert>

// FESystem includes
#include "PostProcess/ElemPostProcessQty.h"

// libMesh includes



ElemPostProcessQty::ElemPostProcessQty(unsigned int el_ID):
  PostProcessQty(),
  calculated_principal_stress(false),
  elem_ID(el_ID)
{

}
	

ElemPostProcessQty::~ElemPostProcessQty()
{


}
  

unsigned int 
ElemPostProcessQty::getElemID() const
{
  return this->elem_ID;
}



void 
ElemPostProcessQty::addStressTensor(TensorValue<double>& value, unsigned int load_case, unsigned int DV_ID)
{
  std::string name;
  name = "Stress";
  
  std::ostringstream stream;
  stream << "_lc_" << load_case;

  if (DV_ID != 0)
    {
      stream <<  "_dv_" << DV_ID;
    }
  
  name += stream.str();

  this->addTensor(name, value);
}



void 
ElemPostProcessQty::addStrainTensor(TensorValue<double>& value, unsigned int load_case, unsigned int DV_ID)
{
  std::string name;
  name = "Strain";
  
  std::ostringstream stream;
  stream << "_lc_" << load_case;

  if (DV_ID != 0)
    {
      stream <<  "_dv_" << DV_ID;
    }
  
  name += stream.str();

  this->addTensor(name, value);
}


void 
ElemPostProcessQty::addMechanicalStrainTensor(TensorValue<double>& value, unsigned int load_case, unsigned int DV_ID)
{
  std::string name;
  name = "MechanicalStrain";
  
  std::ostringstream stream;
  stream << "_lc_" << load_case;

  if (DV_ID != 0)
    {
      stream <<  "_dv_" << DV_ID;
    }
  
  name += stream.str();

  this->addTensor(name, value);
}



TensorValue<double>& 
ElemPostProcessQty::getStressTensor(unsigned int load_case, unsigned int DV_ID)
{
  std::string name;
  name = "Stress";
  
  std::ostringstream stream;
  stream << "_lc_" << load_case;

  if (DV_ID != 0)
    {
      stream <<  "_dv_" << DV_ID;
    }
  
  name += stream.str();

  return this->getTensor(name);
}


TensorValue<double>& 
ElemPostProcessQty::getStrainTensor(unsigned int load_case, unsigned int DV_ID)
{
  std::string name;
  name = "Strain";
  
  std::ostringstream stream;
  stream << "_lc_" << load_case;

  if (DV_ID != 0)
    {
      stream <<  "_dv_" << DV_ID;
    }
  
  name += stream.str();

  return this->getTensor(name);

}

void 
ElemPostProcessQty::calculatePrincipalStressTensor( unsigned int load_case, unsigned int DV_ID)
{
  std::string name;
  name = "PrincipalStress";
  
  std::ostringstream stream;
  stream << "_lc_" << load_case;

  if (DV_ID != 0)
    {
      stream <<  "_dv_" << DV_ID;
    }
  
  name += stream.str();

  // for now, do this for a case without sensitivity
  assert (DV_ID == 0);

  TensorValue<double>& stress = 
    this->getStressTensor(load_case);

  TensorValue<double> value;
  value(0,0) = (stress(0,0) + stress(1,1))*0.5 + 
    pow( pow((stress(0,0) - stress(1,1))*0.5, 2) + stress(0,1)*stress(0,1), 0.5);
  value(1,1) = (stress(0,0) + stress(1,1))*0.5 - 
    pow( pow((stress(0,0) - stress(1,1))*0.5, 2) + stress(0,1)*stress(0,1), 0.5);

  this->addTensor(name, value);

  this->calculated_principal_stress = true;
}



TensorValue<double>& 
ElemPostProcessQty::getPrincipalStressTensor(unsigned int load_case, unsigned int DV_ID)
{
  if (!this->calculated_principal_stress)
    this->calculatePrincipalStressTensor(load_case, DV_ID);

  std::string name;
  name = "PrincipalStress";
  
  std::ostringstream stream;
  stream << "_lc_" << load_case;

  if (DV_ID != 0)
    {
      stream <<  "_dv_" << DV_ID;
    }
  
  name += stream.str();

  return this->getTensor(name);

}



TensorValue<double>& 
ElemPostProcessQty::getMechanicalStrainTensor(unsigned int load_case, unsigned int DV_ID)
{
  std::string name;
  name = "MechanicalStrain";
  
  std::ostringstream stream;
  stream << "_lc_" << load_case;

  if (DV_ID != 0)
    {
      stream <<  "_dv_" << DV_ID;
    }
  
  name += stream.str();

  return this->getTensor(name);
}


double 
ElemPostProcessQty::getStrainEnergyDensity(unsigned int load_case, unsigned int DV_ID)
{
  // get the stress and mechanical strain. 
  
  if (DV_ID == 0)
    {
      TensorValue<double>& stress = this->getStressTensor(load_case);
      TensorValue<double>& strain = this->getMechanicalStrainTensor(load_case);

      // multiply the two to get the strain energy
      double strain_energy_density = 0.0;
      for (unsigned int i=0; i<2; i++)
	for (unsigned int j=0; j<2; j++) 
	  strain_energy_density += stress(i,j)*strain(i,j) ; 
      return strain_energy_density * 0.5;
    }
  else 
    {
      // get the sensitivity and strain
      TensorValue<double>& stress = this->getStressTensor(load_case);
      TensorValue<double>& strain = this->getMechanicalStrainTensor(load_case);
      TensorValue<double>& stress_sens = this->getStressTensor(load_case, DV_ID);
      TensorValue<double>& strain_sens = this->getMechanicalStrainTensor(load_case, DV_ID);

      double strain_energy_density_sens = 0.0;
      for (unsigned int i=0; i<2; i++)
	for (unsigned int j=0; j<2; j++) 
	  strain_energy_density_sens += stress_sens(i,j)*strain(i,j) +
	    stress(i,j) * strain_sens(i,j); 
      
      return strain_energy_density_sens * 0.5;
    }
}





double 
ElemPostProcessQty::getVonMisesStress(unsigned int load_case, unsigned int DV_ID)
{
 
  // for now, this is only for a case without sensitiity. Hence, make sure that
  // DV_ID = 0;
  assert (DV_ID == 0);

  TensorValue<double>& principal = this->getPrincipalStressTensor(load_case);
  
  double von_mises = 0.0;
  von_mises = 
    pow( principal(0,0) - principal(1,1) , 2) +
    principal(0,0)*principal(0,0) + 
    principal(1,1)*principal(1,1);

  von_mises = pow(von_mises*0.5, 0.5);
  
  return von_mises;
}


