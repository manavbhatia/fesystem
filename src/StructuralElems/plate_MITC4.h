// $Id: plate_MITC4.h,v 1.11 2006-10-29 05:09:01 manav Exp $

#ifndef __fesystem_plate_mitc4_h__
#define __fesystem_plate_mitc4_h__

// C++ includes

// FESystem includes
#include "StructuralElems/structural_elem.h"


namespace Discipline
{
  class AnalysisDisciplineBase;
}


#ifndef STRUCTURAL_PLATE_MITC4_QUAD4_ENUM_ID
#define STRUCTURAL_PLATE_MITC4_QUAD4_ENUM_ID 12
#else
#error
#endif


#ifndef STRUCTURAL_PLATE_MITC4_QUAD4_ENUM_NAME
#define STRUCTURAL_PLATE_MITC4_QUAD4_ENUM_NAME "STRUCTURAL_PLATE_MITC4_QUAD4"
#else
#error
#endif


DeclareFESystemElemTypeEnumName(STRUCTURAL_PLATE_MITC4_QUAD4,
                                STRUCTURAL_PLATE_MITC4_QUAD4_ENUM_ID,
                                STRUCTURAL_PLATE_MITC4_QUAD4_ENUM_NAME,
                                QUAD4)




#ifndef PLATE_MITC4_SHEAR_STIFFNESS_FACTOR_ENUM_ID
#define PLATE_MITC4_SHEAR_STIFFNESS_FACTOR_ENUM_ID 1
#else
#error
#endif

#ifndef PLATE_MITC4_SHEAR_STIFFNESS_FACTOR_ENUM_NAME
#define PLATE_MITC4_SHEAR_STIFFNESS_FACTOR_ENUM_NAME "PLATE_MITC4_SHEAR_STIFFNESS_FACTOR"
#else
#error
#endif



namespace FESystemElem
{

  DeclareEnumClass(PlateMITC4QtyEnum);
  
  DeclareEnumName(PLATE_MITC4_SHEAR_STIFFNESS_FACTOR, PlateMITC4QtyEnum,
                  PLATE_MITC4_SHEAR_STIFFNESS_FACTOR_ENUM_ID,
                  PLATE_MITC4_SHEAR_STIFFNESS_FACTOR_ENUM_NAME);
  
  class PlateMITC4: public FESystemElem::StructuralElem
  {
public:
    
    PlateMITC4(Discipline::AnalysisDisciplineBase& discipline);
    
    ~PlateMITC4();
    
    
    /// this function returns a vector of the post process quanties for this element. The quantities will 
    /// consist of the strain, stress and their sensitivities. 
    /// @param vector of load case IDs for which the post process quantities need to be calculated
    /// @param vector of DVs for which the quantities need to be calculated
    virtual std::auto_ptr<ElemPostProcessQty> 
      getElementPostProcessQty(std::vector<unsigned int> , std::vector<DesignData::DesignParameter*>);
    
    
    
    /// this method returns the quantity specified in the arguement list. This has 
    /// been templatized so that the template can be used in derived element types.
    /// This element class will provide its element quantity and calculation routines.
    /// @param name  name of the quantity to be returned
    /// @param qty data structure in which the quantity will be removed
    /// @param side_num number of the side if this is a side quantity
    template <class QtyType >
      void getPlateMITC4Qty(const unsigned int name,
                            QtyType* qty,
                            const unsigned int design_point_enum_ID,
                            const unsigned int domain,
                            const bool sensitivity);
    
protected:
      
      /// inserts the finite element types to be used for this element in the input 
      /// parameters
      virtual void getFETypes(std::vector<FEType>& fetypes);
    
    
    /// inserts the quadrature rules to be used for the different elements
    virtual void getQuadratureRules
      (std::map<FEFamily, std::pair<QuadratureType, Order> >& quadratures);
    
    
    /// this method is an interface to the specific functions to calculate element
    /// quantities. It calls the method that performs the calculations and returns the
    /// quantity in the first arguement. 
    void calculateMITC4ElemQty(DenseMatrix<double>* qty, 
                               const unsigned int name, 
                               const unsigned int design_point_enum_ID,
                               const unsigned int domain);
    
    
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
    
    void calculateShearStiffnessFactor(DenseMatrix<double>* matrix,
                                       const unsigned int design_point_enum_ID);
    
    
    /// this is an  method to calculate the strain operator matrices
    /// @param the vector where all the strain operators will be returned
    /// @param vector of points, defined in the element local coordinate system
    /// where the operators will be calculated
    //virtual void calculateStrainOperator(std::vector<DenseMatrix<double> >* , const std::vector<Point>* );
    
  };
}




template <class QtyType >
void
FESystemElem::PlateMITC4::getPlateMITC4Qty(const unsigned int name,
                                           QtyType* qty,
                                           const unsigned int design_point_enum_ID,
                                           const unsigned int domain,
                                           const bool sensitivity)
{
  // the quantities are not stored in the database.
  assert (qty != NULL);
  qty->zero();
  
  switch (sensitivity)
    {
    case false:      
      this->calculateMITC4ElemQty(qty, name, design_point_enum_ID, domain);
      break;
      
    case true:
    default:
      {
        // this is only for shape sensitivity, since the local quantities specific to this 
        // element are assumed not to be property dependent
        switch (this->sensitivity_parameter)
          {
          case SHAPE_PARAMETER_ENUM_ID:
            {
              // keep going
            }
            break;
            
          default:
            {
              // not for property sensitivity
              Assert(false, ExcInternalError())
            }
          }
        
        static QtyType  base_qty;
        
        // next, obtain the quantity for the perturbed element
        this->calculateMITC4ElemQty(qty, name,
                                      FESystemElem::BASE_PLUS_DELTA_ELEM::num(), domain);
        
        // also, obtain the quantity for the base element. The base quantity also needs to
        // be evaluated at the same load case, dof_value for sensitivity calculation
        base_qty.zero();
        this->calculateMITC4ElemQty(&base_qty, name, FESystemElem::BASE_ELEM::num(), domain);
        
        qty->add(-1.0, base_qty);
        qty->scale(1.0/this->perturbation);
      }
      break;
    }
}






#endif // __fesystem_plate_mitc4_h__
