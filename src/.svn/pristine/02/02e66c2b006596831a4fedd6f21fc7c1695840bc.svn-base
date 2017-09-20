// $Id: ForcedConvection1D.h,v 1.3 2007-01-05 03:08:10 manav Exp $

#ifndef __fesystem_forced_convection_1d_h__
#define __fesystem_forced_convection_1d_h__

// C++ includes


// FESystem includes
#include "thermal_elem.h"

// libMesh includes


namespace Discipline
{
  class AnalysisDisciplineBase;
}


#ifndef THERMAL_FORCED_CONVECTION_EDGE2_ENUM_ID
#define THERMAL_FORCED_CONVECTION_EDGE2_ENUM_ID 19
#else
#error
#endif


#ifndef THERMAL_FORCED_CONVECTION_EDGE2_ENUM_NAME
#define THERMAL_FORCED_CONVECTION_EDGE2_ENUM_NAME "THERMAL_FORCED_CONVECTION_EDGE2"
#else
#error
#endif


DeclareFESystemElemTypeEnumName(THERMAL_FORCED_CONVECTION_EDGE2,
                                THERMAL_FORCED_CONVECTION_EDGE2_ENUM_ID,
                                THERMAL_FORCED_CONVECTION_EDGE2_ENUM_NAME,
                                QUAD4)


namespace FESystemElem
{
  
  class ForcedConvection_1D: public ThermalElem
  {
public:
    ForcedConvection_1D(Discipline::AnalysisDisciplineBase& discipline);
    
    ~ForcedConvection_1D();
    
    /// @returns the number of dofs in the element
    virtual unsigned int getNDofs();

protected:
      
      /**
      *	function to calculate the element capacitance matrix. The first arguement is the matrix
       *	that will be resized to the appropriate size, and filled with the matrix elements. If 
       *	sensitivity of the matrix wrt to a property is to be calculated, then second arguement 
       *	with the name of the property should be provided. 
       */
      virtual void calculate_C(DenseMatrix<double>* matrix, 
                               const unsigned int design_point_enum_ID,
                               bool sensitivity_calculation);
    
    /**
      *	function to calculate the element conductance matrix. The first arguement is the matrix
     *	that will be resized to the appropriate size, and filled with the matrix elements. If 
     *	sensitivity of the matrix wrt to a property is to be calculated, then second arguement 
     *	with the name of the property should be provided. 
     */
    virtual void calculate_K_c(DenseMatrix<double>* matrix, 
                               const unsigned int design_point_enum_ID,
                               bool sensitivity_calculation);
    
    /**
      *	function to calculate the conductance of element conductance matrix. 
     * The first arguement is the matrix
     *	that will be resized to the appropriate size, and filled with the matrix elements. 
     * This quantity exists only when the material properties are temperature dependent. 
     * Otherwise, it is zero
     */
    virtual void calculate_K_c_Jac(DenseMatrix<double>* matrix,
                                   const unsigned int design_point_enum_ID);
    
    
    /**
      *	function to calculate the element capacitance matrix contribution due to 
     *	surface convection load. The first arguement is the matrix
     *	that will be resized to the appropriate size, and filled with the matrix elements. 
     *  The second arguement should be the vector of surface convection loads for the element.
     *	If sensitivity of the matrix wrt to a property is to be calculated, then third arguement 
     *	with the name of the property should be provided. 
     */
    virtual void calculate_K_h(DenseMatrix<double>* matrix, 
                               const unsigned int design_point_enum_ID,
                               bool sensitivity_calculation);
    
    /**
      *	function to calculate the element surface  heat flux load vector.
     *	The first arguement is the vector
     *	that will be resized to the appropriate size, and filled with the vector elements. 
     *  The second arguement should be the vector of surface flux loads for the element.
     *	If sensitivity of the vector wrt to a property is to be calculated, then third arguement 
     *	with the name of the property should be provided. 
     */
    virtual void calculate_F_qsurf(DenseVector<double>* vector, 
                                   const unsigned int design_point_enum_ID,
                                   bool sensitivity_calculation);
    
    /**
      *	function to calculate the element volume heat load vector.
     *	 The first arguement is the vector
     *	that will be resized to the appropriate size, and filled with the vector elements. 
     *  The second arguement should be the vector of volume heat loads for the element.
     *	If sensitivity of the vector wrt to a property is to be calculated, then third arguement 
     *	with the name of the property should be provided. 
     */
    virtual void calculate_F_Qvol(DenseVector<double>* vector, 
                                  const unsigned int design_point_enum_ID,
                                  bool sensitivity_calculation);
    
    
    /**
      *	function to calculate the element surface convection load vector.
     *  The first arguement is the vector
     *	that will be resized to the appropriate size, and filled with the vector elements. 
     *  The second arguement should be the vector of surface convection loads for the element.
     *	If sensitivity of the vector wrt to a property is to be calculated, then third arguement 
     *	with the name of the property should be provided. 
     */
    virtual void calculate_F_h(DenseVector<double>* vector, 
                               const unsigned int design_point_enum_ID,
                               bool sensitivity_calculation);
    
    
    /**
      *	function to calculate the element emitted radiation load vector.
     *	The first arguement is the vector
     *	that will be resized to the appropriate size, and filled with the vector elements. 
     *  The third arguement should be the vector of emitted radiation load vector for the element.
     *	If sensitivity of the vector wrt to a property is to be calculated, then third arguement 
     *	with the name of the property should be provided. 
     */
    virtual void calculate_F_sigma(DenseVector<double>* vector, 
                                   const unsigned int design_point_enum_ID,
                                   bool sensitivity_calculation);
    
    /// calculated and returns the jacobian contribution from the radiation loads
    virtual void calculate_F_sigma_Jac(DenseMatrix<double>* matrix,
                                       const unsigned int design_point_enum_ID);
    
    /// upwind weighting factor
    double upwind_factor;
    
  };
  
}


inline
unsigned int FESystemElem::ForcedConvection_1D::getNDofs()
{
  return this->getNNodes();
}

#endif //__fesystem_forced_convection_1d_h__
