// $Id: FEInterpolationAnalysis.h,v 1.3 2006-09-05 20:41:48 manav Exp $

#ifndef __fesystem_interpolation_analysis_h__
#define __fesystem_interpolation_analysis_h__

// C++ includes 


// FESystem includes 
#include "Discipline/AnalysisDisciplineBase.h"

// libMesh includes
#include "numerics/sparse_matrix.h"
#include "numerics/numeric_vector.h"
#include "geom/elem.h"

// Forward Decleration
namespace FESystem
{
  class FESystemController;
}


namespace FESystemElem
{
  class FESystemElemBase;
}

#ifndef FE_INTERPOLATION_DISCIPLINE_ENUM_ID
#define FE_INTERPOLATION_DISCIPLINE_ENUM_ID 3
#else
#error
#endif

#ifndef FE_INTERPOLATION_DISCIPLINE_ENUM_NAME
#define FE_INTERPOLATION_DISCIPLINE_ENUM_NAME "FE_INTERPOLATION_DISCIPLINE"
#else
#error
#endif


DeclareEnumName(FE_INTERPOLATION_DISCIPLINE, Discipline::AnalysisDisciplineEnum,
                FE_INTERPOLATION_DISCIPLINE_ENUM_ID,
                FE_INTERPOLATION_DISCIPLINE_ENUM_NAME);


/// This class provides methods to calculate the matrices and vectors related to the 
/// thermal discipline. 
/// this is a problem with this class. It adds a variable to the mesh/dofmap, in addition
/// to the ones already added by the structural/thermal disciplines, hence distorts the
/// dof numbering, etc. Hence, these variables should be added to an independent dofmap,
/// without messing with the thermal/structural mesh. It basically needs a separate 
/// dofmap structure.
class FEInterpolationAnalysis: public Discipline::AnalysisDisciplineBase
{
 public:
	
  /// constructor. It takes two arguements,
  /// the arguement is a pointer to the AnalysisDriver that owns this object
  FEInterpolationAnalysis(FESystem::FESystemController& controller,
                          Discipline::InterpolationDisciplineInfo& info);
	
  /// destructor
  ~FEInterpolationAnalysis();
	
  /// returns if same system matrix can be used for all load cases
  virtual inline bool sameSystemMatrixForLoadCases();
  
  /// function will initialize the matrices and vectors for this anlysis
  void initMatricesAndVectors();

	
 protected:
	
    /// function to add analysis variables
    virtual void addVariablesForAnalysis();
  
  
  /// function will return a string name for the quantity
  /// @param quantity enum name
  /// @param bool saying if the quantity is for sensitiivty or not
  std::string getStringNameForQty(const unsigned int qty_enum_ID , bool sensitivity);
  
  
  /// function will fill the input vector with the element quantities required for the calculation
  /// of the specified global quantity. This is an abstract function and will be implemented in the
  /// derived class
  /// @param global quantitiy
  /// @param element pointer
  /// @param vector in which the quantities will be added
  void getElemQty(const unsigned int qty_enum_ID ,
                  FESystemElem::FESystemElemBase* elem,
                  std::vector<DenseMatrix<double>* >*	qty_vector);
  
  
  /// function will fill the input vector with the element quantities required for 
  /// the calculation
  /// of the specified global quantity. This is an abstract function and will be implemented 
  /// in the
  /// derived class
  /// @param global quantitiy
  /// @param element pointer
  /// @param vector in which the quantities will be added
  void getElemQty(const unsigned int qty_enum_ID ,
                  FESystemElem::FESystemElemBase* elem,
                  std::vector<DenseVector<double>* >*	qty_vector);
  
  
  /// function to calculate the conductance matrix. First arguement is the matrix to store 
  ///  the result in, and the second arguement is a bool to tell the function if sensitivity
  /// is desired
  void calculateK(SparseMatrix<double>& , bool = false);
  
  /// function to calculate the Jacobian matrix. First arguement is the matrix to store 
  ///  the result in
  void calculateJac(SparseMatrix<double>& );
  
  /// function to calculate the force vector. First arguement is the vector to store 
  ///  the result in, and the second arguement is a bool to tell the function if sensitivity
  /// is desired
  void calculateF(NumericVector<double>& , bool = false);
};



inline 
bool
FEInterpolationAnalysis::sameSystemMatrixForLoadCases()
{
  switch (this->analysis_type_enum_ID)
    {
    case LINEAR_ANALYSIS_ENUM_ID:
      return true;
      break;
      
    case NONLINEAR_ANALYSIS_ENUM_ID:
      abort();
      break;
      
    default:
      abort();
      break;
    }
}



#endif // __fesystem_interpolation_analysis_h__
