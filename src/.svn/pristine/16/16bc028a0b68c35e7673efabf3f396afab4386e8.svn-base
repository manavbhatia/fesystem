// $Id: AnalysisDisciplineBase.h,v 1.14.6.6 2008-06-03 05:19:32 manav Exp $

#ifndef __fesystem_analysis_discipline_base_h__
#define __fesystem_analysis_discipline_base_h__

// C++ includes 
#include <memory>
#include <vector>
#include <set>

// FESystem includes
#include "Utilities/NameEnumHandler.h"
#include "FESystem/FESystemController.h"
#include "Mesh/MeshList.h"
#include "Mesh/FEMeshData.h"
#include "FESystem/FESystemElem.h"
#include "AnalysisDriver/AnalysisDriver.h"
#include "Utilities/LoadedSolution.h"
#include "Utilities/ParallelUtility.h"

// libMesh includes
#include "base/libmesh_common.h"
#include "fe/fe_type.h"
#include "numerics/sparse_matrix.h"
#include "numerics/numeric_vector.h"
#include "numerics/petsc_matrix.h"

// Forward declerations
class ElementDataStorage;

namespace FESystemElem
{
  class FESystemElemBase;
}

namespace Discipline
{
  class DisciplineInfo;
}

namespace Driver
{
  class AnalysisDriver;
}

namespace FESystemDatabase
{
  class DataInfoBase;
}

namespace MeshDS
{
  class FEMeshData;
}

namespace Solution
{
  class SolutionBase;
}

namespace FESystemUtility
{
  class LoadedSolution;
}


namespace Loads
{
  class NodalLoadCombination;
  class VolumeLoadCombination;
  class SurfaceLoadCombination;
  class DirichletBoundaryConditionDataInfo;
}

class FEType;
class Elem;
template <typename T> class DenseMatrix;
template <typename T> class DenseVector;
template <typename T> class SparseMatrix;
template <typename T> class NumericVector;



#ifndef LINEAR_ANALYSIS_ENUM_ID 
#define LINEAR_ANALYSIS_ENUM_ID 1
#else
#error
#endif

#ifndef LINEAR_ANALYSIS_ENUM_NAME 
#define LINEAR_ANALYSIS_ENUM_NAME "LINEAR_ANALYSIS"
#else
#error
#endif

#ifndef NONLINEAR_ANALYSIS_ENUM_ID 
#define NONLINEAR_ANALYSIS_ENUM_ID 2
#else
#error
#endif

#ifndef NONLINEAR_ANALYSIS_ENUM_NAME 
#define NONLINEAR_ANALYSIS_ENUM_NAME "NONLINEAR_ANALYSIS"
#else
#error
#endif



#ifndef STEADY_STATE_SOLUTION_ENUM_ID 
#define STEADY_STATE_SOLUTION_ENUM_ID 1
#else
#error
#endif

#ifndef STEADY_STATE_SOLUTION_ENUM_NAME 
#define STEADY_STATE_SOLUTION_ENUM_NAME "STEADY_STATE_SOLUTION"
#else
#error
#endif


#ifndef TRANSIENT_SOLUTION_ENUM_ID 
#define TRANSIENT_SOLUTION_ENUM_ID 2
#else
#error
#endif

#ifndef TRANSIENT_SOLUTION_ENUM_NAME 
#define TRANSIENT_SOLUTION_ENUM_NAME "TRANSIENT_SOLUTION"
#else
#error
#endif



namespace Discipline
{
  
  DeclareEnumClass(AnalysisDisciplineEnum);
  DeclareEnumClass(DisciplineAnalysisTypeEnum);

  DeclareEnumName(LINEAR_ANALYSIS, DisciplineAnalysisTypeEnum,
                  LINEAR_ANALYSIS_ENUM_ID,
                  LINEAR_ANALYSIS_ENUM_NAME);
  
  DeclareEnumName(NONLINEAR_ANALYSIS, DisciplineAnalysisTypeEnum,
                  NONLINEAR_ANALYSIS_ENUM_ID,
                  NONLINEAR_ANALYSIS_ENUM_NAME);
  
  
  DeclareEnumClass(SolutionTransientNatureTypeEnum);

  DeclareEnumName(STEADY_STATE_SOLUTION, SolutionTransientNatureTypeEnum,
                  STEADY_STATE_SOLUTION_ENUM_ID,
                  STEADY_STATE_SOLUTION_ENUM_NAME);

  DeclareEnumName(TRANSIENT_SOLUTION, SolutionTransientNatureTypeEnum,
                  TRANSIENT_SOLUTION_ENUM_ID,
                  TRANSIENT_SOLUTION_ENUM_NAME);

  /// This class provides methods to calculate the matrices and vectors related to the 
  /// thermal discipline. 
  class AnalysisDisciplineBase
    {
public:
      
      /// constructor. It takes two arguements, first is a reference to the FESystemController 
      /// that contains the analysis data like the meshlist, propertylist, load-database, etc. 
      /// and the second is a pointer to the AnalysisDisciplineDriver that owns this object
      AnalysisDisciplineBase(FESystem::FESystemController& controller,
                             const Discipline::DisciplineInfo& info);
      
      /// destructor
      virtual ~AnalysisDisciplineBase();
      
      /// function to clear all data structures
      void clear();

      /// @returns the disipline enum ID
      inline unsigned int getDisciplineEnumID() const;
      
      /// @returns the discipline enum name
      inline std::string getDisciplineEnumName() const;
      
      /// sets the solution kind
      inline void attachSolution(Solution::SolutionBase* sol);
      
      /// @returns a reference to the solution
      const Solution::SolutionBase& getSolution() const;
      
      /// @returns the discipline info for this object
      inline const Discipline::DisciplineInfo& getDisciplineInfo() const;
      
      /// attaches an analysis driver to this discipline
      void attachAnalysisDriver(Driver::AnalysisDriver *driver);
      
      /// @returns a FESystemDatabase::DataInfo pointer, which defines the quantity 
      /// passed in the parameter, for the current analysis
      virtual std::auto_ptr<FESystemDatabase::DataInfoBase> getDataInfoForQty
      (const unsigned int qty_enum_ID) = 0;

      /// @returns the type of analysis for this discipline, i.e. if it is linear or nonlinear
      unsigned int getAnalysisTypeEnumID() const;

      /// @returns the type of analysis for this discipline, i.e. if it is linear or nonlinear
      const std::string& getAnalysisTypeEnumName() const;
      

      inline FESystem::FESystemController& getFESystemController();

      inline Driver::AnalysisDriver& getAnalysisDriver();
      
      inline unsigned int getAnalysisMeshID() const;
      
      /// returns the Mesh for this analysis 
      inline const MeshDS::FEMesh& getAnalysisMesh() const;
      
      /// returns the MeshData for this analysis
      inline const MeshDS::FEMeshData& getAnalysisMeshData() const;
      
      /// returns the DOFMap object associated with the mesh
      inline const MeshDS::FEDofMap& getAnalysisDofMap() const;
      
      /// returns the dummy elem list
      inline const std::set<Elem*>& getAnalysisDummyElemSet() const;
      
      /// @returns the number of variables in the analysis
      inline unsigned int getNVars() const;
      
      /// @returns the number of degrees of freedom in the analysis
      inline unsigned int getNDofs() const;

      /// @returns the eigenproblem kind for the discipline. This has to be 
      /// implemented for each derived class.
      virtual unsigned int getEigenProblemKindEnumID() = 0;

      
      /// @returns the name of variable \p i.
      inline const std::string & getVariableName(const unsigned int i) const;
      
      
      /// @returns the variable number assoicated with
      /// the user-specified variable named \p var.
      inline unsigned short int getVariableNumber (const std::string& var) const;
      
      
      /// @returns the finite element type variable number \p i.
      inline const FEType & getVariableType (const unsigned int i) const;
      
      
      /// @returns the finite element type for variable \p var.
      inline const FEType & getVariableType (const std::string& var) const;
      
      
      /// checks if the analysis is dependent on temperature
      inline bool checkPropertyDependenceOnTemperature() const;
      
      /// returns the enum ID of the boundary condition for this disicipline
      virtual std::auto_ptr<Loads::DirichletBoundaryConditionDataInfo> 
        getBoundaryConditionLoadInfo() const = 0;
      

      /// @returns the order of the transient nature of the problem
      virtual unsigned int getTransientSystemOrder() = 0;
      
      /// @returns true if the specified matrix has exists for the discipline
      virtual bool disciplineHasMatrix(const unsigned int qty_enum_ID) const = 0;
      
      
      /// this function fills in the matrix/vector with the quantity referred to by
      /// the second arguement
      void fillQuantity(std::map<unsigned int, double>& real_map, 
                        std::map<unsigned int, SparseMatrix<double>* >& matrix_map, 
                        std::map<unsigned int, NumericVector<double>* >& vector_map);
          
      
      /// this function will iterate over all the elements, calculate the post process
      /// quantity, and 
      /// store them in the post process quantity database 
      void calculatePostProcessQty(const std::vector<unsigned int>& load_cases);
      
      /// @returns the perturbed elem for the given elem ID and design param
      Elem* getPerturbedElemForShapeParameter(const unsigned int elem_ID, 
                                              DesignData::DesignParameter* param);
        
      
      /// this method returns a boolean saying if the same system matrix can be used 
      /// for each load case
      virtual bool sameSystemMatrixForLoadCases() = 0;
      
      /// this will initialize the matrix to the sparsity pattern of the discipline's mesh
      virtual void initMatrix(SparseMatrix<double>& mat);


      /// this will initialize the vector to the sparsity pattern of the discipline's mesh
      virtual void initVector(NumericVector<double>& vec);
      
      /// this function will return a vector containing the dof values for the element 
      /// nodes in a var major format.
      /// @param load case ID
      void getElemDofValues(const Elem* elem,
                            DenseVector<double>& qty, 
                            const unsigned int transient_order,
                            const bool sensitivity);
      
      /// @returns the number of dofs on the discipline
      unsigned int nDofs() const;
      
      
      /// @returns the number of dofs on the current processor
      unsigned int nLocalDofs() const;
      
      
protected:
        
      
      /// Adds the variable \p var to the list of variables in this analysis
      void addVariable(const std::string& var, const FEType& type);
      
      /// set the number of systems for each DOFObject
      void setNSystemsForDOFObjects();
      
      
      /// function to add analysis variables. This is an abstract function and
      /// will be implemented in the derived class
      virtual void addAnalysisVariables() = 0;

      /// adds variables and distributes the dofs
      void addAnalysisVariablesAndInitialize();
      
      
      /// function creates elements for the specific discipline. This is an abstract 
      /// function and will be implemented in the derived class. The input arguements
      /// are 
      /// @param element ID
      /// @param fesys_elem pointer to the fesystem element that will be used for analysis
      inline void getFESystemElem(const unsigned int elem_kind_enum_ID, 
                                  FESystemElem::FESystemElemBase*& fesys_elem); 
      
      
      /// function will fill the input vector with the element quantities required for the calculation
      /// of the specified global quantity. This is an abstract function and will be implemented in the
      /// derived class
      /// @param global quantitiy
      /// @param element pointer
      /// @param vector in which the quantities will be added
      virtual void getElemQty(FESystemElem::FESystemElemBase* elem,
                              std::map<unsigned int, double>& real_elem_data,
                              std::map<unsigned int, DenseMatrix<double> >& mat_elem_data,
                              std::map<unsigned int, DenseVector<double> >& vec_elem_data,
                              const bool sensitivity,
                              const unsigned int DV_ID) = 0;
      
      
      ///  this function iterates over all elements and calculates the global quantity
      /// @param global data
      /// @param quantity to be calculates
      /// @param whether sensitivity is required
      void calculateGlobalQty
      (std::map<unsigned int, double >& global_real_data_map,
       std::map<unsigned int, SparseMatrix<double>* >& global_matrix_data_map,
       std::map<unsigned int, NumericVector<double>* >& global_vector_data_map,
       bool sensitivity);
      
      
      /// this method resizes the matrices and vectors if needed.
      void resizeElemMatricesForElem
        ( FESystemElem::FESystemElemBase* fesys_elem,
         std::map<unsigned int, DenseMatrix<double> >& elem_matrix_data_map,
         std::map<unsigned int, DenseVector<double> >& elem_vector_data_map);
        
      /// function will add the element data to the appropriate locations in the 
      /// global matrix.
      /// @param global matrix
      /// @param libmesh elem pointer
      /// @param vector containing element quantities
      inline void addElemQtyToGlobalData
        (const Elem* elem, 
         std::map<unsigned int, double>& global_real_data_map,
         std::map<unsigned int, double>& elem_real_data_map,
         std::map<unsigned int, SparseMatrix<double>*>& global_matrix_data_map,
         std::map<unsigned int, DenseMatrix<double> >& elem_matrix_data_map,
         std::map<unsigned int, NumericVector<double>*>& global_vector_data_map,
         std::map<unsigned int, DenseVector<double> >& elem_vector_data_map);
      
      
      /// this gets the loads from the load database
      void getElemLoads
	(const unsigned int elem_ID,  const Elem* elem,
	 const unsigned int load_case, 
	 const bool transient, const double time, 
	 const bool sensitivity, const unsigned int DV_ID,
	 const std::map<unsigned int, unsigned int>& volume_loads,
	 const std::map<unsigned int, unsigned int>& surface_loads,
	 const std::map<unsigned int, std::pair<unsigned int, unsigned int> >& nodal_loads,
	 std::map<unsigned int, const Loads::VolumeLoadCombination*>& volume_load_map,
	 std::map<unsigned int, std::map<unsigned int, const Loads::NodalLoadCombination*> >& nodal_load_map,
	 std::map<unsigned int, std::map<unsigned int, const Loads::SurfaceLoadCombination*> >& surface_load_map);
      
      
      /// this method fills the vector with the solution vector necessary for filling the quantity
      /// for the current analysis. It checks the current status of the discipline, and 
      /// apropriately loads the solution vector
      virtual const NumericVector<double>& getSolutionVector(const unsigned int transient_order,
                                                             const bool sensitivity) = 0;
          
      /// function to calculate the conductance matrix. First arguement is the matrix to store 
      ///  the result in, and the second arguement is a bool to tell the function if sensitivity
      /// is desired
      virtual bool calculateTransientMatrix(const unsigned int order, bool = false) = 0;

      /// function to calculate the conductance matrix. First arguement is the matrix to store 
      ///  the result in, and the second arguement is a bool to tell the function if sensitivity
      /// is desired
      virtual bool calculateK(bool = false) = 0;
      
      /// function to calculate the Jacobian matrix. First arguement is the matrix to store 
      ///  the result in, and the second arguement is a bool to tell the function if sensitivity
      /// is desired
      virtual bool calculateJac( ) = 0;
      
      /// function to calculate the A matrix for eigen problem analysis
      virtual bool calculateEigenProblemAMatrix(bool sensitivity = false) = 0;
      
      /// function to calculate the B matrix for eigen problem analysis
      virtual bool calculateEigenProblemBMatrix(bool sensitivity = false) = 0;
      
      /// function to calculate the force vector. First arguement is the vector to store 
      ///  the result in, and the second arguement is a bool to tell the function if sensitivity
      /// is desired
      virtual bool calculateF( bool = false) = 0;
      
      /// this method will perform some discipline specific operations once the 
      /// global quantities have been calculated
      virtual void performAdditionalOperationsOnGlobalQuantities
        (std::map<unsigned int, SparseMatrix<double>* >& matrix_map,
         std::map<unsigned int, NumericVector<double>* >& vector_map) = 0;
      
      FESystem::FESystemController& fesystem_controller;
      
      const Discipline::DisciplineInfo& discipline_info;

      /// dependence on temperature
      bool temperature_dependent_property;
      
      /// discipline name enum ID
      const unsigned int discipline_enum_ID;
      
      /// type of the analysis, i.e. linear, nonlinear
      const unsigned int analysis_type_enum_ID;
      
      /// solution kind enum ID
      Solution::SolutionBase* solution_base;
      
      /// pointer to the analysis driver for this discipline
      Driver::AnalysisDriver *analysis_driver;
      
      /// ID of mesh
      const unsigned int mesh_ID;
      
      /// mesh data structure
      MeshDS::FEMesh* mesh;
      
      /// mesh data data-structure
      MeshDS::FEMeshData* mesh_data;
      
      /// dof map data-structure
      MeshDS::FEDofMap* dof_map;
      
      /// set of elements that do not need to be processed for element quantity 
      /// calculations
      std::set<Elem*>* dummy_element_set;

      /// this is the vector that is loaded from the database to provide the dofs
      /// to the elements when needed from a previously calculated solution
      FESystemUtility::LoadedSolution loaded_solution;
      
      /// The names of the variables associated with this analysis
      std::vector<System::Variable*> _variables;
            
      /// map of fesystem elements for this discipline. These are initialized using the
      /// geometric elements from libMesh to calculate element quantities. The key used here
      /// is the element type enum ID.
      typedef std::map<unsigned int, FESystemElem::FESystemElemBase*> FESystemElemMap;
      FESystemElemMap fesystem_elem_map;
    };
  
  
  
}


inline
FESystem::FESystemController& 
Discipline::AnalysisDisciplineBase::getFESystemController() 
{
  return this->fesystem_controller;
}



inline
Driver::AnalysisDriver& 
Discipline::AnalysisDisciplineBase::getAnalysisDriver() 
{
  Assert(this->analysis_driver != NULL, ExcEmptyObject());
  return *(this->analysis_driver);
}


inline
void
Discipline::AnalysisDisciplineBase::attachSolution(Solution::SolutionBase* sol)
{
  Assert(sol != NULL, ExcEmptyObject());
  
  this->solution_base = sol;
}


inline 
unsigned int
Discipline::AnalysisDisciplineBase::getDisciplineEnumID() const
{
  return this->discipline_enum_ID;
}



inline 
std::string
Discipline::AnalysisDisciplineBase::getDisciplineEnumName() const
{
  return Discipline::AnalysisDisciplineEnum::enumName(this->discipline_enum_ID);
}



inline
const Discipline::DisciplineInfo& 
Discipline::AnalysisDisciplineBase::getDisciplineInfo() const
{
  return this->discipline_info;
}



inline
const MeshDS::FEDofMap& 
Discipline::AnalysisDisciplineBase::getAnalysisDofMap() const
{
  return *(this->dof_map);
}


inline
unsigned int
Discipline::AnalysisDisciplineBase::getAnalysisMeshID() const
{
  return this->mesh_ID;
}




inline
const MeshDS::FEMesh& 
Discipline::AnalysisDisciplineBase::getAnalysisMesh() const
{
  return *(this->mesh);
}



inline
const MeshDS::FEMeshData&
Discipline::AnalysisDisciplineBase::getAnalysisMeshData() const
{
  return *(this->mesh_data);
}



inline 
const std::set<Elem*>& 
Discipline::AnalysisDisciplineBase::getAnalysisDummyElemSet() const
{
  return *(this->dummy_element_set);
}



inline
unsigned int 
Discipline::AnalysisDisciplineBase::getNVars() const
{
  return this->_variables.size();
}



inline
unsigned int
Discipline::AnalysisDisciplineBase::getNDofs() const
{
  return this->dof_map->n_dofs();
}





inline 
const std::string& 
Discipline::AnalysisDisciplineBase::getVariableName(const unsigned int i) const
{
  assert (i < this->getNVars());
	
  return this->_variables[i]->name();
}




inline
unsigned short int 
Discipline::AnalysisDisciplineBase::getVariableNumber(const std::string& var) const
{
  for (unsigned int i=0; i<this->_variables.size(); i++)
    if (this->_variables[i]->name() == var)
      return this->_variables[i]->number();
  
  Assert(false, ExcInternalError());
  return NULL;
}


inline
const FEType& 
Discipline::AnalysisDisciplineBase::getVariableType(const unsigned int i) const
{
  assert (i < this->getNVars());

  return this->_variables[i]->type();
}



inline
const FEType& 
Discipline::AnalysisDisciplineBase::getVariableType(const std::string& var) const
{
  for (unsigned int i=0; i<this->_variables.size(); i++)
    if (this->_variables[i]->name() == var)
      return this->_variables[i]->type();

  Assert(false, ExcInternalError());
  return FEType();
}



inline
bool
Discipline::AnalysisDisciplineBase::checkPropertyDependenceOnTemperature() const
{
  return this->temperature_dependent_property;
}



inline void 
Discipline::AnalysisDisciplineBase::getFESystemElem(const unsigned int elem_kind_enum_ID,
                                                    FESystemElem::FESystemElemBase*& fesys_elem)
{
  
  // check if this elem kind exists in the map
  unsigned int elem_exists = this->fesystem_elem_map.count(elem_kind_enum_ID);
  
  // if it does not exist, create and insert the elem
  if (elem_exists == 0)
    {
    FESystemElem::FESystemElemBase* local_fesys_elem = 
    FESystemElem::createFESystemElem(elem_kind_enum_ID, *this).release();
    
    bool insert_success = 
      this->fesystem_elem_map.insert(Discipline::AnalysisDisciplineBase::FESystemElemMap::
                                     value_type(elem_kind_enum_ID, local_fesys_elem)).second;
     
    // this should never be false, but still, a check is introduced
    Assert(insert_success, ExcInternalError());
    }
  
  // now return the element
  fesys_elem = this->fesystem_elem_map[elem_kind_enum_ID];
}



inline
void 
Discipline::AnalysisDisciplineBase::addElemQtyToGlobalData
(const Elem* elem, 
 std::map<unsigned int, double>& global_real_data_map,
 std::map<unsigned int, double>& elem_real_data_map,
 std::map<unsigned int, SparseMatrix<double>*>& global_matrix_data_map,
 std::map<unsigned int, DenseMatrix<double> >& elem_matrix_data_map,
 std::map<unsigned int, NumericVector<double>*>& global_vector_data_map,
 std::map<unsigned int, DenseVector<double> >& elem_vector_data_map)
{
  Assert(global_matrix_data_map.size() == elem_matrix_data_map.size(), 
         ExcInternalError());
  Assert(global_vector_data_map.size() == elem_vector_data_map.size(),
         ExcInternalError());

  // get the dof_indices
  static std::vector<unsigned int> dof_indices;
  
  this->dof_map->dof_indices(elem, dof_indices);
  
  static std::map<unsigned int, double>::iterator real_it, real_end;
  static std::map<unsigned int, SparseMatrix<double>*>::iterator mat_it, mat_end;
  static std::map<unsigned int, NumericVector<double>*>::iterator vec_it, vec_end;

  real_it = global_real_data_map.begin();
  real_end = global_real_data_map.end();
  
  mat_it = global_matrix_data_map.begin();
  mat_end = global_matrix_data_map.end();

  vec_it = global_vector_data_map.begin();
  vec_end = global_vector_data_map.end();

  for ( ; real_it != real_end; real_it++)
    real_it->second += elem_real_data_map[real_it->first];
  
  for ( ; mat_it != mat_end; mat_it++)
    mat_it->second->add_matrix(elem_matrix_data_map[mat_it->first], dof_indices);
  
  
  for ( ; vec_it != vec_end; vec_it++)
    vec_it->second->add_vector(elem_vector_data_map[vec_it->first], dof_indices);
}






inline
void
Discipline::AnalysisDisciplineBase::getElemDofValues(const Elem* elem,
                                                     DenseVector<double>& qty,
                                                     const unsigned int transient_order,
                                                     const bool sensitivity)
{
  static unsigned int n_nodes, dof, n_vars;
  n_nodes = elem->n_nodes();
  n_vars = this->getNVars();
  
  // create a quantity and resize it
  if (qty.size() != (n_nodes * n_vars))
    qty.resize(n_nodes * n_vars);
  
  const NumericVector<double>& sol = this->getSolutionVector(transient_order, sensitivity);
  
  // all computations inside the element are in a variable major mode
  for (unsigned int var_it =0; var_it < n_vars; var_it++)
    for (unsigned int node_it =0; node_it < n_nodes; node_it++)
      {
      dof = elem->get_node(node_it)->dof_number(0,var_it,0);
      qty(var_it * n_nodes + node_it) = sol(dof);
      }
}


#endif // __fesystem_analysis_discipline_base_h__
