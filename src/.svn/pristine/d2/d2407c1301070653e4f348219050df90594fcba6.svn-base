// $Id: FEInterpolationElem.h,v 1.5 2006-10-23 23:42:35 manav Exp $

#ifndef __fesystem_interpolation_elem_h__
#define __fesystem_interpolation_elem_h__

// C++ includes
#include <memory>

// FESystem includes
#include "FESystem/FESystemElem.h"
#include "AnalysisDriver/AnalysisDriver.h"


// libMesh includes
#include "numerics/dense_matrix.h"
#include "numerics/dense_vector.h"


// Forward decleration
class Load;
class FEInterpolationAnalysis;


#ifndef FE_INTERPOLATION_ELEM_K_MATRIX_ENUM_ID
#define FE_INTERPOLATION_ELEM_K_MATRIX_ENUM_ID 1
#else
#error
#endif

#ifndef FE_INTERPOLATION_ELEM_K_MATRIX_ENUM_NAME
#define FE_INTERPOLATION_ELEM_K_MATRIX_ENUM_NAME "FE_INTERPOLATION_ELEM_K_MATRIX"
#else
#error
#endif


#ifndef FE_INTERPOLATION_ELEM_F_VECTOR_ENUM_ID
#define FE_INTERPOLATION_ELEM_F_VECTOR_ENUM_ID 2
#else
#error
#endif

#ifndef FE_INTERPOLATION_ELEM_F_VECTOR_ENUM_NAME
#define FE_INTERPOLATION_ELEM_F_VECTOR_ENUM_NAME "FE_INTERPOLATION_ELEM_F_VECTOR"
#else
#error
#endif

DeclareEnumClass(FEInterpolationElemQtyEnum);

DeclareEnumName(FE_INTERPOLATION_ELEM_K_MATRIX, FEInterpolationElemQtyEnum,
                FE_INTERPOLATION_ELEM_K_MATRIX_ENUM_ID, FE_INTERPOLATION_ELEM_K_MATRIX_ENUM_NAME);

DeclareEnumName(FE_INTERPOLATION_ELEM_F_VECTOR, FEInterpolationElemQtyEnum,
                FE_INTERPOLATION_ELEM_F_VECTOR_ENUM_ID, FE_INTERPOLATION_ELEM_F_VECTOR_ENUM_NAME);


/**
*	a base class defining the heat conduction elements
 */

class FEInterpolationElem: public FESystemElem::FESystemElemBase
{
public:
  FEInterpolationElem(const unsigned int dim, 
                      const unsigned int elem_enum_ID,
                      Discipline::AnalysisDisciplineBase& discipline);
	
  virtual ~FEInterpolationElem();
  
	
  /// function will get the matrix or vector quantity for this element
  /// specified by the name in the arguement. This function will obtain the
  /// required information for calculating the quantity. Hence, the user 
  /// need not worry about passsing them as arguements. The function will 
  /// check for the current load case from the analysis driver and get the 
  /// loads if needed. For nonlinear loads, it will procure the displacements
  /// by itself. 
  void  getElementAssembledQty(const unsigned int quantity,
                               DenseMatrix<double>* return_qty);
	
  
  /// function will get the matrix or vector quantity for this element
  /// specified by the name in the arguement. This function will obtain the
  /// required information for calculating the quantity. Hence, the user 
  /// need not worry about passsing them as arguements. The function will 
  /// check for the current load case from the analysis driver and get the 
  /// loads if needed. For nonlinear loads, it will procure the displacements
  /// by itself. 
  void  getElementAssembledQty(const unsigned int quantity,
                               DenseVector<double>* return_qty);
  
  
  /// this function will return the sensitivity quantity by procuring the 
  /// required information by itself. It will check for the current load case and design
  /// variable and perform the operations by itself.
  void  getElementAssembledQtySensitivity(const unsigned int quantity,
                                          DenseMatrix<double>* return_qty);
  
  
  /// this function will return the sensitivity quantity by procuring the 
  /// required information by itself. It will check for the current load case and design
  /// variable and perform the operations by itself.
  void  getElementAssembledQtySensitivity(const unsigned int quantity,
                                          DenseVector<double>* return_qty);
  
  
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
  void calculateAssembledQty(DenseMatrix<double>* ,
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
  
  
//  /// returns a string name for the quantity. If it is a side quantity, the last arguement 
//  /// must contain the side number
//  std::string getStringNameForFEInterpolationElemQty(const unsigned int qty,
//                                                     const unsigned int domain);
//	
//	
//  /// returns a string name for the quantity sensitivity. The second arguement is the 
//  /// DV number, and the last arguement is the optional side number, needed only for
//  /// side quantity
//  std::string getStringNameForFEInterpolationElemQtySensitivity(const unsigned int qty, 
//                                                                const unsigned int DV,
//                                                                const unsigned int domain);
	
  
  /// calculates the interpolation stiffness matrix
  virtual void calculate_K(DenseMatrix<double>* matrix, 
                           const unsigned int design_point_enum_ID,
                           bool sensitivity_calculation);
	
	
  /// calculates the interpolation force vector
  virtual void calculate_F(DenseVector<double>* vector, 
                           const unsigned int design_point_enum_ID,
                           bool sensitivity_calculation);
  
  
};


#endif // __interpolation_elem_h__
