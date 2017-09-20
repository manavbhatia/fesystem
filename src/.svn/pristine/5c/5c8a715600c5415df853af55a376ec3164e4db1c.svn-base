// $Id: membrane.h,v 1.12 2006-10-23 23:44:01 manav Exp $

#ifndef __fesystem_membrane_h__
#define __fesystem_membrane_h__

// C++ includes

// FESystem includes
#include "StructuralElems/structural_elem.h"


namespace Discipline
{
  class AnalysisDisciplineBase;
}



namespace FESystemElem
{
  /**
  *	a class definig the membrane structural element
   */
  
  class Membrane: public StructuralElem
  {
public:
    Membrane(const unsigned int dim, 
             const unsigned int elem_enum_ID,
             Discipline::AnalysisDisciplineBase& discipline);
    
    ~Membrane();
    
    
    /// this function returns a vector of the post process quanties for this element. The quantities will 
    /// consist of the strain, stress and their sensitivities. 
    /// @param vector of load case IDs for which the post process quantities need to be calculated
    /// @param vector of DVs for which the quantities need to be calculated
    virtual std::auto_ptr<ElemPostProcessQty> 
      getElementPostProcessQty(std::vector<unsigned int> , 
                               std::vector<DesignData::DesignParameter*>);
    
    
protected:
      
      /// inserts the finite element types to be used for this element in the input 
      /// parameters
      virtual void getFETypes(std::vector<FEType>& fetypes);
    
    
    /// inserts the quadrature rules to be used for the different elements
    virtual void getQuadratureRules
      (std::map<FEFamily, std::pair<QuadratureType, Order> >& quadratures);
    
      virtual void calculate_M(DenseMatrix<double>*, 
                               const unsigned int design_point_enum_ID,
                               bool sensitivity_calculation);
	   
    virtual void calculate_K(DenseMatrix<double>*, 
			     const unsigned int design_point_enum_ID,
                             bool sensitivity_calculation);
	   
    /// calculates and returns the stiffness matrix
    virtual void calculate_K_G(DenseMatrix<double>*, 
                               const unsigned int design_point_enum_ID,
			       bool sensitivity_calculation);

    virtual void calculate_F_T(DenseVector<double>*, 
                               const unsigned int design_point_enum_ID,
                               bool sensitivity_calculation);
	   
	   
    /// this is an  method to calculate the strain operator matrices
    /// @param the vector where all the strain operators will be returned
    /// @param vector of points, defined in the element local coordinate system where the operators will be 
    /// calculated
    //virtual void calculateStrainOperator(std::vector<DenseMatrix<double> >* , const std::vector<Point>* );
    
    
  };
}

#endif // __fesystem_membrane_h__
