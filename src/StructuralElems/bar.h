// $Id: bar.h,v 1.20 2006-10-29 05:09:01 manav Exp $

#ifndef __fesystem_bar_h__
#define __fesystem_bar_h__

// C++ includes

// FESystem includes
#include "StructuralElems/structural_elem.h"


namespace Discipline
{
  class AnalysisDisciplineBase;
}



#ifndef STRUCTURAL_BAR_EDGE2_ENUM_ID
#define STRUCTURAL_BAR_EDGE2_ENUM_ID 6
#else
#error
#endif


#ifndef STRUCTURAL_BAR_EDGE2_ENUM_NAME
#define STRUCTURAL_BAR_EDGE2_ENUM_NAME "STRUCTURAL_BAR_EDGE2"
#else
#error
#endif


DeclareFESystemElemTypeEnumName(STRUCTURAL_BAR_EDGE2,
                                STRUCTURAL_BAR_EDGE2_ENUM_ID,
                                STRUCTURAL_BAR_EDGE2_ENUM_NAME,
                                EDGE2)



namespace FESystemElem
{
  
  /**
  *	a class definig the spring structural element
   */
  
  class Bar: public FESystemElem::StructuralElem
  {
public:
    Bar(Discipline::AnalysisDisciplineBase& discipline);
    
    ~Bar();
    
    /// this function returns a vector of the post process quanties for this element. 
    /// The quantities will 
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
	   
	   
    //	virtual bool calculate_F_B(DenseVector<double>*, 
    //							   Property::PropertyName = Property::INVALID_PROPERTY_NAME);
    
    /// this is an  method to calculate the strain operator matrices
    /// @param the vector where all the strain operators will be returned
    /// @param vector of points, defined in the element local coordinate system where the operators will be 
    /// calculated
    //virtual void calculateStrainOperator(std::vector<DenseMatrix<double> >* , const std::vector<Point>* );
  };
}

#endif // __fesystem_bar_h__
