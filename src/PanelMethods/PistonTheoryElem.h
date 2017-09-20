/*
 *  PistonTheoryElem.h
 *  fesystem_xcode
 *
 *  Created by Manav Bhatia on 9/21/08.
 *  Copyright 2008 Manav Bhatia. All rights reserved.
 *
 */

// $Id:$

#ifndef __fesystem_piston_theory_elem_h__
#define __fesystem_piston_theory_elem_h__

// C++ includes


// FESystem includes
#include "FESystem/FESystemElem.h"
#include "Utilities/NameEnumHandler.h"
#include "Properties/ElemDataCard.h"
#include "Discipline/AnalysisDisciplineBase.h"
#include "Mesh/FEMeshData.h"

// libMesh includes
#include "numerics/dense_matrix.h"
#include "numerics/dense_vector.h"


// Forward declerations
class Load;
namespace Discipline
{
  class AnalysisDisciplineBase;
}


#ifndef PISTON_THEORY_K_MATRIX_ENUM_ID
#define PISTON_THEORY_K_MATRIX_ENUM_ID 1
#else
#error
#endif

#ifndef PISTON_THEORY_K_MATRIX_ENUM_NAME
#define PISTON_THEORY_K_MATRIX_ENUM_NAME "PISTON_THEORY_K_MATRIX"
#else
#error
#endif


#ifndef PISTON_THEORY_C1_MATRIX_ENUM_ID
#define PISTON_THEORY_C1_MATRIX_ENUM_ID 2
#else
#error
#endif

#ifndef PISTON_THEORY_C1_MATRIX_ENUM_NAME
#define PISTON_THEORY_C1_MATRIX_ENUM_NAME "PISTON_THEORY_C2_MATRIX"
#else
#error
#endif



namespace FESystemElem {
  
  DeclareEnumClass(PistonTheoryElemQtyEnum);
  
  DeclareEnumName(PISTON_THEORY_K_MATRIX, FESystemElem::PistonTheoryElemQtyEnum,
                  PISTON_THEORY_K_MATRIX_ENUM_ID, PISTON_THEORY_K_MATRIX_ENUM_NAME);

  DeclareEnumName(PISTON_THEORY_C1_MATRIX, FESystemElem::PistonTheoryElemQtyEnum,
                  PISTON_THEORY_C1_MATRIX_ENUM_ID, PISTON_THEORY_C1_MATRIX_ENUM_NAME);

  
  class PistonTheoryElem : public FESystemElemBase {
    
  public:
    /// constructor
    PistonTheoryElem(const unsigned int dim, 
                     const unsigned int elem_enum_ID,
                     Discipline::AnalysisDisciplineBase& discipline);
    
    /// destructor
    virtual ~PistonTheoryElem();
    
    
    /// @returns the number of dofs for this element
    virtual unsigned int getNDofs();
    
    /// @returns the transient system order for the elem equations
    virtual unsigned int getTransientSystemOrder() const;
    
    /// function will get the matrix or vector quantity for this element
    /// specified by the name in the arguement. 
    void getElementAssembledQty(const unsigned int name_enum_ID,
                                DenseVector<double>* qty);
    
    
    /// function will get the matrix or vector quantity for this element
    /// specified by the name in the arguement. 
    void getElementAssembledQty(const unsigned int name_enum_ID,
                                DenseMatrix<double>* qty);
    
    
    /// this function will return the sensitivity quantity by procuring the 
    /// required information by itself. It will check for the current load case and design
    /// variable and perform the operations by itself.
    void  getElementAssembledQtySensitivity(const unsigned int name_enum_ID,
                                            DenseVector<double>* qty);
    
    
    /// this function will return the sensitivity quantity by procuring the 
    /// required information by itself. It will check for the current load case and design
    /// variable and perform the operations by itself.
    void  getElementAssembledQtySensitivity(const unsigned int name_enum_ID,
                                            DenseMatrix<double>* qty);
    
    
    
    /// sets the vector of the fluid flow in global axis
    void setFluidFlowVector(const Point& vector);
    
    /// clears the data structures of this element
    virtual void clearLoadAndDofs();
    
    /// this function returns the element post process quantity
    /// @param vector of load case IDs for which the quantitis should be calculated
    /// @param vector of DVs for which the quantities should be calculated
    virtual std::auto_ptr<ElemPostProcessQty>
    getElementPostProcessQty(std::vector<unsigned int> load_cases,
                             std::vector<DesignData::DesignParameter*> DV_vector);
    

  protected:
    
    /// inserts the finite element types to be used for this element in the input 
    /// parameters
    virtual void getFETypes(std::vector<FEType>& fetypes);
    
    
    /// inserts the quadrature rules to be used for the different elements
    virtual void getQuadratureRules
    (std::map<FEFamily, std::pair<QuadratureType, Order> >& quadratures);

    /// function to direct the call to appropriate calculation routine based on the 
    /// quantity seeked. The arguements are:
    /// @param  pointer of quantity to return the data in
    /// @param quantity name
    /// @param bool specifying if sensitivity has been asked
    void calculateAssembledQty(DenseMatrix<double>* qty,
                               const unsigned int qty_enum_ID,
                               const unsigned int design_point_enum_ID,
                               bool sensitivity_calculation);
    
    
    /// function to direct the call to appropriate calculation routine based on the 
    /// quantity seeked. The arguements are:
    /// @param  pointer of quantity to return the data in
    /// @param quantity name
    /// @param bool specifying if sensitivity has been asked
    void calculateAssembledQty(DenseVector<double>* qty,
                               const unsigned int qty_enum_ID,
                               const unsigned int design_point_enum_ID,
                               bool sensitivity_calculation);
    
    
    /// function will return the nodal temperature vector in the first arguement of the function
    /// @param vector in which the temperatures will be returned
    /// @param bool indicating if sensitivity is being sought
    /// @param DV ID
    void extractNodalTemperatureVectorFromLoads(DenseVector<double>& temp_vec, 
                                                bool sensitivity);
    
    
    /// calculates the transformation matrix for the structural element
    void getTransformationMatrix(DenseMatrix<double>* matrix,
                                 const unsigned int design_point_enum_ID,
                                 bool shape_sensitivity);
    
    /**
     *	function will calculate the matrix of the element,
     *	and fill the matrix passed as the arguement. The second arguement is the 
     *	name of the property, which, if passed to this function, will make the 
     *	function calculate the sensitivity of the qty wrt the property
     */
    virtual void calculate_C(DenseMatrix<double>*, 
                             const unsigned int design_point_enum_ID,
                             bool sensitivity_calculation);
    
    /// calculates and returns the element stiffness matrix
    virtual void calculate_K(DenseMatrix<double>*, 
                             const unsigned int design_point_enum_ID,
                             bool sensitivity_calculation);
    
    /// This is a dummy function, and does not do anything
    virtual void initPropertyCard();
    
    /// This is a dummy function and does not do anything
    virtual void clearPropertyCardInitialization();
    
    /// the vector of fluid flow 
    DenseVector<double> fluid_flow_vector;
          
  };
}



inline
unsigned int FESystemElem::PistonTheoryElem::getNDofs()
{
  return 6 * (this->getNNodes());
}


inline
unsigned int FESystemElem::PistonTheoryElem::getTransientSystemOrder() const
{
  return 1;
}


#endif // __fesystem_piston_theory_elem_h__

