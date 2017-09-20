// $Id: StructuralAnalysis.C,v 1.19.6.10 2008-08-21 00:37:18 manav Exp $

// C++ includes
#include <sstream>

// FESystem includes
#include "Discipline/StructuralAnalysis.h"
#include "Discipline/StructuralDisciplineInfo.h"
#include "FESystem/FESystemController.h"
#include "FESystem/AnalysisCase.h"
#include "Mesh/MeshList.h"
#include "Mesh/FEMeshData.h"
#include "Database/GlobalDataStorage.h"
#include "Loads/LoadDatabase.h"
#include "AnalysisDriver/AnalysisDriver.h"
#include "StructuralElems/structural_elem.h"
#include "AnalysisDriver/LinearAnalysisDriver.h"
#include "AnalysisDriver/EigenProblemAnalysisDriver.h"
#include "AnalysisDriver/TransientAnalysisDriver.h"
#include "Loads/LoadDatabase.h"
#include "Loads/LoadDataInfo.h"
#include "Loads/LoadCombination.h"
#include "FESystem/FESystemExceptions.h"
#include "Solvers/EigenSolver.h"
#include "Solutions/LinearStressSolution.h"
#include "Solutions/StructuralVibrationEigenSolution.h"
#include "Solutions/LinearizedBucklingEigenSolution.h"
#include "Solutions/TransientStructuralSolution.h"
#include "Solutions/AeroelasticitySolution.h"
#include "Solutions/NonlinearStaticStructuralSolution.h"


// libMesh includes
#include "numerics/sparse_matrix.h"
#include "numerics/numeric_vector.h"
#include "numerics/dense_matrix.h"
#include "numerics/dense_vector.h"


Discipline::StructuralAnalysis::StructuralAnalysis
(FESystem::FESystemController& controller,
 const Discipline::StructuralDisciplineInfo& info):
Discipline::AnalysisDisciplineBase(controller, info)
{
  this->addAnalysisVariablesAndInitialize();
}



Discipline::StructuralAnalysis::~StructuralAnalysis()
{
	
}



void
Discipline::StructuralAnalysis::addAnalysisVariables()
{  
  FEType fetype(FIRST, LAGRANGE);
  // add a temperature variable to the analysis driver
  this->addVariable("u", fetype);
  this->addVariable("v", fetype);
  this->addVariable("w", fetype);
  this->addVariable("theta_x", fetype);
  this->addVariable("theta_y", fetype);
  this->addVariable("theta_z", fetype);
}





unsigned int
Discipline::StructuralAnalysis::getEigenProblemKindEnumID()
{
  unsigned int solution_kind = this->solution_base->getSolutionEnumID();
  switch(solution_kind)
  {
    case STRUCTURAL_VIBRATION_EIGEN_SOLUTION_ENUM_ID:
    case LINEARIZED_BUCKLING_EIGEN_SOLUTION_ENUM_ID:
    {
      return Solver::GENERALIZED_HERMITIAN_EIGENPROBLEM::num();
    }
      break;
      
    default:
      Assert(false, FESystemExceptions::ExcEnumCaseNotImplemented
             (Solution::SolutionEnum::enumName(solution_kind)));
  }
  
  // shold not get here
  return FESystemNumbers::InvalidID;
}


std::auto_ptr<Loads::DirichletBoundaryConditionDataInfo> 
Discipline::StructuralAnalysis::getBoundaryConditionLoadInfo() const
{
  Assert(this->analysis_driver != NULL, ExcInternalError());
  
  std::auto_ptr<Loads::DirichletBoundaryConditionDataInfo> return_ptr;
  
  return_ptr.reset(new Loads::DirichletBoundaryConditionDataInfo);
  return_ptr->setLoadCaseID(this->analysis_driver->getCurrentLoadCase());
  return_ptr->setLoadNameEnumID(DISPLACEMENT_BOUNDARY_CONDITION::num());
  return_ptr->setLoadClassEnumID(DIRICHLET_BOUNDARY_CONDITION::num());
  
  // check if the analysis is time dependent
  switch (this->solution_base->getDisciplineTransientNatureEnumID
          (Discipline::STRUCTURAL_DISCIPLINE::num()))
  {
    case STEADY_STATE_SOLUTION_ENUM_ID:
      // nothing to be done here
      break;
      
    case TRANSIENT_SOLUTION_ENUM_ID:
    {
      double time_val = dynamic_cast<Driver::TransientAnalysisDriver*>
      (this->analysis_driver)->getCurrentAnalysisTime();
      return_ptr->setTime(time_val);
    }
      break;
      
    default:
      // at this point of time, no other analysis needs the specification of the elem
      // dof vector
      Assert(false, ExcInternalError());
      
  }
  
  return return_ptr;
}





bool 
Discipline::StructuralAnalysis::calculateTransientMatrix(const unsigned int order,
                                                         bool sensitivity)
{
  // this is to avoid getting warnings about unused parameters
  (void) sensitivity;
  (void) order;
  // for now, everything will be calculated. the below logic does not take into account 
  // nonlinear analysis. that will be done later.
  return true;
  
  //   unsigned int load_case = this->analysis_driver->getCurrentLoadCase();
	
  //   // get the name by which the quantity will be referred
  //   std::string name = this->getStringNameForQty(Driver::SYSTEM_MATRIX::num(), sensitivity);
  
  //   // check if the K matrix exists in the database
  //   bool K_matrix_present = 
  //     this->fesystem_controller.global_data_storage->checkIfMatrixExists
  //     (load_case, name);
  
  //   // if K is not present, then it needs to be recalculated
  //   return !K_matrix_present;
}


bool 
Discipline::StructuralAnalysis::calculateK(bool sensitivity)
{
  // this is to avoid getting warnings about unused parameters
  (void) sensitivity;
  // for now, everything will be calculated. the below logic does not take into account 
  // nonlinear analysis. that will be done later.
  return true;
  
  //   unsigned int load_case = this->analysis_driver->getCurrentLoadCase();
	
  //   // get the name by which the quantity will be referred
  //   std::string name = this->getStringNameForQty(Driver::SYSTEM_MATRIX::num(), sensitivity);
  
  //   // check if the K matrix exists in the database
  //   bool K_matrix_present = 
  //     this->fesystem_controller.global_data_storage->checkIfMatrixExists
  //     (load_case, name);
  
  //   // if K is not present, then it needs to be recalculated
  //   return !K_matrix_present;
}





bool 
Discipline::StructuralAnalysis::calculateJac()
{
  // there still is no nonlinear structural analysis capability yet.
  Assert(false, ExcInternalError());
  
  // returning a bogus value to avoid warnings of no return
  return false;
}




bool
Discipline::StructuralAnalysis::calculateEigenProblemAMatrix(bool sensitivity) 
{
  // this is to avoid getting warnings about unused parameters
  (void) sensitivity;
  // for now, everything will be calculated. the below logic does not take into account 
  // nonlinear analysis. that will be done later.
  return true;
  
  //   unsigned int load_case = this->analysis_driver->getCurrentLoadCase();
	
  //   // get the name by which the quantity will be referred
  //   std::string name = this->getStringNameForQty
  //     (Driver::EIGENPROBLEM_A_MATRIX::num(), sensitivity);
	
  //   // check if the K matrix exists in the database
  //   bool matrix_present = 
  //     this->fesystem_controller.global_data_storage->checkIfMatrixExists
  //     (load_case, name);
	
  //   return !matrix_present;
}




bool
Discipline::StructuralAnalysis::calculateEigenProblemBMatrix(bool sensitivity)
{
  // this is to avoid getting warnings about unused parameters
  (void) sensitivity;
  // for now, everything will be calculated. the below logic does not take into account 
  // nonlinear analysis. that will be done later.
  return true;
  
  //   unsigned int load_case = this->analysis_driver->getCurrentLoadCase();
	
  //   // get the name by which the quantity will be referred
  //   std::string name = this->getStringNameForQty
  //     (Driver::EIGENPROBLEM_B_MATRIX::num(), sensitivity);
	
  //   // check if the K matrix exists in the database
  //   bool matrix_present = 
  //     this->fesystem_controller.global_data_storage->checkIfMatrixExists
  //     (load_case, name);
  
  //   return !matrix_present;
}



bool 
Discipline::StructuralAnalysis::calculateF(bool sensitivity)
{
  // this is to avoid getting warnings about unused parameters
  (void) sensitivity;
  // for now, everything will be calculated. the below logic does not take into account 
  // nonlinear analysis. that will be done later.
  return true;
  
  //   unsigned int load_case = this->analysis_driver->getCurrentLoadCase();
	
  //   // get the name by which the quantity will be referred
  //   std::string name = this->getStringNameForQty(Driver::FORCE_VECTOR::num(), sensitivity);
	
  //   // check if the vector exists in the database
  //   bool Force_vector_present = 
  //     this->fesystem_controller.global_data_storage->checkIfVectorExists
  //     (load_case, name);
	
  //   return !Force_vector_present;
}




void 
Discipline::StructuralAnalysis::getElemQty
(FESystemElem::FESystemElemBase* elem,
 std::map<unsigned int, double>& real_elem_data,
 std::map<unsigned int, DenseMatrix<double> >& mat_elem_data,
 std::map<unsigned int, DenseVector<double> >& vec_elem_data,
 const bool sensitivity,
 const unsigned int DV_ID)
{
  // make sure that valid pointer have been given 
  assert (elem != NULL);
	
  static DenseMatrix<double> matrix;
  static DenseVector<double> vector, elem_dofs, elem_dof_sensitivity;
  static std::vector<DenseVector<double>*> elem_dof_vec(3), elem_dof_sens_vec(3);
  
  // get loads for the current load case
  static unsigned int load_case, elem_ID;
  // the loads will be stored in this map
  static std::map<unsigned int, const Loads::VolumeLoadCombination*> volume_loads;
  static std::map<unsigned int, std::map<unsigned int, const Loads::SurfaceLoadCombination*> > surface_loads;
  static std::map<unsigned int, std::map<unsigned int, const Loads::NodalLoadCombination*> > nodal_loads;
  
  static std::map<unsigned int, const Loads::VolumeLoadCombination*> volume_load_sensitivity;
  static std::map<unsigned int, std::map<unsigned int, const Loads::SurfaceLoadCombination*> > 
  surface_load_sensitivity;
  static std::map<unsigned int, std::map<unsigned int, const Loads::NodalLoadCombination*> > nodal_load_sensitivity;
  
  static std::map<unsigned int, unsigned int> volume_load_map;
  static std::map<unsigned int, unsigned int> surface_load_map;
  static std::map<unsigned int, std::pair<unsigned int, unsigned int> > nodal_load_map;
  
  static bool created_load_maps = false;
  
  if (!created_load_maps)
    {
      nodal_load_map[NODAL_TEMPERATURE::num()] =  std::pair<unsigned int, unsigned int>(NODAL_POINT_LOAD::num(), 1);
      surface_load_map[SURFACE_PRESSURE::num()] = SCALAR_SURFACE_LOAD::num();
//      volume_load_map[SELF_WEIGHT::num()] = VECTOR_VOLUME_LOAD::num();
      created_load_maps = true;
    }
  
  elem_ID = elem->getID();
  load_case = this->analysis_driver->getCurrentLoadCase();
  
  // set the transient nature of this discipline
  static double time_val = 0.0;
  static bool transient_nature = false;
  
  if (this->solution_base->getDisciplineTransientNatureEnumID(Discipline::STRUCTURAL_DISCIPLINE::num()) 
      == Discipline::TRANSIENT_SOLUTION::num())
    {
      transient_nature = true;
      time_val = dynamic_cast<Driver::TransientAnalysisDriver*>(this->analysis_driver)->getCurrentAnalysisTime();
    }
  
  
  this->getElemLoads(elem_ID, elem->getElem(FESystemElem::BASE_ELEM::num()), 
                     load_case, transient_nature, time_val,
                     false, FESystemNumbers::InvalidID,
                     volume_load_map, surface_load_map, nodal_load_map,
                     volume_loads, nodal_loads, surface_loads);
  
  if (!sensitivity)
    elem->setElementLoads(&volume_loads, &surface_loads, &nodal_loads,
                          NULL, NULL, NULL);
  else
    {
      this->getElemLoads(elem_ID, elem->getElem(FESystemElem::BASE_ELEM::num()), 
                         load_case, transient_nature, time_val,
                         true, DV_ID,
                         volume_load_map, surface_load_map, nodal_load_map,
                         volume_load_sensitivity, nodal_load_sensitivity, surface_load_sensitivity);
      elem->setElementLoads(&volume_loads, &surface_loads, &nodal_loads,
                            &volume_load_sensitivity, &surface_load_sensitivity,
                            &nodal_load_sensitivity);
    }
  
  // iterate over the matrix quantities and get the data
  static std::map<unsigned int, DenseMatrix<double> >::iterator mat_it, mat_end;
  mat_it = mat_elem_data.begin();
  mat_end = mat_elem_data.end();
  
  for (; mat_it != mat_end; mat_it++)
    {
      mat_it->second.zero();
      
      switch(mat_it->first)
      {
        case TRANSIENT_C2_MATRIX_ENUM_ID:
        {
          elem->getElementAssembledQty(FESystemElem::STRUCTURAL_M_MATRIX::num(),
                                       &matrix);
          mat_it->second += matrix;
        }
          break;
          
        case SYSTEM_MATRIX_ENUM_ID:
        {
          elem->getElementAssembledQty(FESystemElem::STRUCTURAL_K_MATRIX::num(),
                                       &matrix);
          mat_it->second += matrix;

          switch(this->solution_base->getSolutionEnumID())
          {
            case AEROELASTICITY_SOLUTION_ENUM_ID:
            {
              if (dynamic_cast<Solution::AeroelasticitySolution*>
                  (this->solution_base)->ifIncludeGeometricEffects())
                {
                  elem_dofs.zero();
                  this->getElemDofValues(elem->getElem(FESystemElem::BASE_ELEM::num()),
                                         elem_dofs, 0, false);
                  
                  elem_dof_vec[0] = &elem_dofs;
                  elem_dof_vec[1] = NULL;
                  elem_dof_vec[2] = NULL;
                  
                  elem_dof_sens_vec[0] = NULL;
                  elem_dof_sens_vec[1] = NULL;
                  elem_dof_sens_vec[2] = NULL;
                  
                  elem->setDofValues(elem_dof_vec, elem_dof_sens_vec);
                  
                  matrix.zero();
                  
                  elem->getElementAssembledQty(FESystemElem::STRUCTURAL_K_G_MATRIX::num(),
                                               &matrix);
                  mat_it->second += matrix;
                  elem->clearDofValues();
                }
              
            }
              break;
          }
        }
          break;
          
          
        case JACOBIAN_MATRIX_ENUM_ID:
        {
          Assert(false, ExcInternalError());
        }
          break;
          

        case SYSTEM_MATRIX_SENSITIVITY_ENUM_ID:
        {
          // get K sensitivity
          elem->getElementAssembledQtySensitivity(FESystemElem::STRUCTURAL_K_MATRIX::num(),
                                                  &matrix);
          mat_it->second += matrix;
          
          switch(this->solution_base->getSolutionEnumID())
          {
            case AEROELASTICITY_SOLUTION_ENUM_ID:
            {
              Assert(false, ExcInternalError());
              
              if (dynamic_cast<Solution::AeroelasticitySolution*>
                  (this->solution_base)->ifIncludeGeometricEffects())
                {
                  elem_dofs.zero();
                  this->getElemDofValues(elem->getElem(FESystemElem::BASE_ELEM::num()),
                                         elem_dofs, 0, false);
                  
                  elem_dof_vec[0] = &elem_dofs;
                  elem_dof_vec[1] = NULL;
                  elem_dof_vec[2] = NULL;
                  
                  elem_dof_sens_vec[0] = NULL;
                  elem_dof_sens_vec[1] = NULL;
                  elem_dof_sens_vec[2] = NULL;
                  
                  elem->setDofValues(elem_dof_vec, elem_dof_sens_vec);
                  
                  matrix.zero();
                  
                  elem->getElementAssembledQty(FESystemElem::STRUCTURAL_K_G_MATRIX::num(),
                                               &matrix);
                  mat_it->second += matrix;
                  elem->clearDofValues();
                }
              
            }
              break;
          }
          
        }
          break;
          
        case EIGENPROBLEM_A_MATRIX_ENUM_ID:
        {
          switch(this->solution_base->getSolutionEnumID())
          {
            case STRUCTURAL_VIBRATION_EIGEN_SOLUTION_ENUM_ID:
            {
              elem->getElementAssembledQty(FESystemElem::STRUCTURAL_K_MATRIX::num(),
                                           &matrix);
              mat_it->second += matrix;
              
              if (dynamic_cast<Solution::StructuralVibrationEigenSolution*>
                  (this->solution_base)->ifIncludeGeometricEffects())
                {
                  elem_dofs.zero();
                  this->getElemDofValues(elem->getElem(FESystemElem::BASE_ELEM::num()),
                                         elem_dofs, 0, false);
                  
                  elem_dof_vec[0] = &elem_dofs;
                  elem_dof_vec[1] = NULL;
                  elem_dof_vec[2] = NULL;
                  
                  elem_dof_sens_vec[0] = NULL;
                  elem_dof_sens_vec[1] = NULL;
                  elem_dof_sens_vec[2] = NULL;
                  
                  elem->setDofValues(elem_dof_vec, elem_dof_sens_vec);
                 
                  matrix.zero();
                  
                  elem->getElementAssembledQty(FESystemElem::STRUCTURAL_K_G_MATRIX::num(),
                                               &matrix);
                  mat_it->second += matrix;
                  elem->clearDofValues();
                }
            }
              break;
              
            case LINEARIZED_BUCKLING_EIGEN_SOLUTION_ENUM_ID:
            {
              elem->getElementAssembledQty(FESystemElem::STRUCTURAL_K_MATRIX::num(),
                                           &matrix);
              mat_it->second += matrix;
            }
              break;
              
            default:
              abort();
          }
        }
          break;
          
          
        case EIGENPROBLEM_B_MATRIX_ENUM_ID:
        {
          switch(this->solution_base->getSolutionEnumID())
          {
            case STRUCTURAL_VIBRATION_EIGEN_SOLUTION_ENUM_ID:
            {
              elem->getElementAssembledQty(FESystemElem::STRUCTURAL_M_MATRIX::num(),
                                           &matrix);
              mat_it->second += matrix;
            }
              break;
              
            case LINEARIZED_BUCKLING_EIGEN_SOLUTION_ENUM_ID:
            {
              elem_dofs.zero();
              this->getElemDofValues(elem->getElem(FESystemElem::BASE_ELEM::num()),
                                     elem_dofs, 0, false);
              
              elem_dof_vec[0] = &elem_dofs;
              elem_dof_vec[1] = NULL;
              elem_dof_vec[2] = NULL;
              
              elem_dof_sens_vec[0] = NULL;
              elem_dof_sens_vec[1] = NULL;
              elem_dof_sens_vec[2] = NULL;
              
              elem->setDofValues(elem_dof_vec, elem_dof_sens_vec);
              
              elem->getElementAssembledQty(FESystemElem::STRUCTURAL_K_G_MATRIX::num(),
                                           &matrix);
              matrix.scale(-1.0);
              mat_it->second += matrix;
              elem->clearDofValues();
            }
              break;
              
            default:
              abort();
          }
        }
          break;
          
        case EIGENPROBLEM_A_MATRIX_SENSITIVITY_ENUM_ID:
        {
          switch(this->solution_base->getSolutionEnumID())
          {
            case STRUCTURAL_VIBRATION_EIGEN_SOLUTION_ENUM_ID:
            {
              elem->getElementAssembledQtySensitivity(FESystemElem::STRUCTURAL_K_MATRIX::num(),
                                                      &matrix);
              mat_it->second += matrix;
              
              if (dynamic_cast<Solution::StructuralVibrationEigenSolution*>
                  (this->solution_base)->ifIncludeGeometricEffects())
                {
                  elem_dofs.zero();
                  elem_dof_sensitivity.zero();
                  this->getElemDofValues(elem->getElem(FESystemElem::BASE_ELEM::num()),
                                         elem_dofs, 0, false);
                  this->getElemDofValues(elem->getElem(FESystemElem::BASE_ELEM::num()),
                                         elem_dof_sensitivity, 0, true);
                  
                  
                  elem_dof_vec[0] = &elem_dofs;
                  elem_dof_vec[1] = NULL;
                  elem_dof_vec[2] = NULL;
                  
                  elem_dof_sens_vec[0] = &elem_dof_sensitivity;
                  elem_dof_sens_vec[1] = NULL;
                  elem_dof_sens_vec[2] = NULL;
                  
                  elem->setDofValues(elem_dof_vec, elem_dof_sens_vec);
                  
                  matrix.zero();
                  
                  elem->getElementAssembledQtySensitivity
                  (FESystemElem::STRUCTURAL_K_G_MATRIX::num(), &matrix);
                  mat_it->second += matrix;
                  
                  elem->clearDofValues();
                }
            }
              break;
              
            case LINEARIZED_BUCKLING_EIGEN_SOLUTION_ENUM_ID:
            {
              elem->getElementAssembledQtySensitivity(FESystemElem::STRUCTURAL_K_MATRIX::num(),
                                                      &matrix);
              mat_it->second += matrix;
            }
              break;
              
            default:
              abort();
          }
        }
          break;
          
          
        case EIGENPROBLEM_B_MATRIX_SENSITIVITY_ENUM_ID:
        {
          switch(this->solution_base->getSolutionEnumID())
          {
            case STRUCTURAL_VIBRATION_EIGEN_SOLUTION_ENUM_ID:
            {
              elem->getElementAssembledQtySensitivity(FESystemElem::STRUCTURAL_M_MATRIX::num(),
                                                      &matrix);
              mat_it->second += matrix;
            }
              break;
              
            case LINEARIZED_BUCKLING_EIGEN_SOLUTION_ENUM_ID:
            {
              elem_dofs.zero();
              elem_dof_sensitivity.zero();
              this->getElemDofValues(elem->getElem(FESystemElem::BASE_ELEM::num()),
                                     elem_dofs, 0, false);
              this->getElemDofValues(elem->getElem(FESystemElem::BASE_ELEM::num()),
                                     elem_dof_sensitivity, 0, true);
              
              elem_dof_vec[0] = &elem_dofs;
              elem_dof_vec[1] = NULL;
              elem_dof_vec[2] = NULL;
              
              elem_dof_sens_vec[0] = &elem_dof_sensitivity;
              elem_dof_sens_vec[1] = NULL;
              elem_dof_sens_vec[2] = NULL;
              
              elem->setDofValues(elem_dof_vec, elem_dof_sens_vec);
                            
              elem->getElementAssembledQtySensitivity
              (FESystemElem::STRUCTURAL_K_G_MATRIX::num(), &matrix);
              matrix.scale(-1.0);
              mat_it->second += matrix;
              
              elem->clearDofValues();
            }
              break;
              
            default:
              abort();
          }
        }
          break;
          
        default:
          abort();
          break;
      }
    }
  
  
  
  // iterate over the matrix quantities and get the data
  static std::map<unsigned int, DenseVector<double> >::iterator vec_it, vec_end;
  vec_it = vec_elem_data.begin();
  vec_end = vec_elem_data.end();
  
  for (; vec_it != vec_end; vec_it++)
    {
      vec_it->second.zero();
      
      switch(vec_it->first)
      {
        case FORCE_VECTOR_ENUM_ID:
        {
          if (nodal_loads.count(NODAL_TEMPERATURE::num()) > 0)
            {
              elem->getElementAssembledQty(FESystemElem::STRUCTURAL_F_T_VECTOR::num(), &vector);
              vec_it->second += vector;
            }
          
          if (surface_loads.count(SURFACE_PRESSURE::num()) > 0)
            {
              elem->getElementAssembledQty(FESystemElem::STRUCTURAL_F_PRESSURE_VECTOR::num(), 
                                           &vector);
              vec_it->second += vector;
            }
          
//          if (volume_loads.count(SELF_WEIGHT::num()) > 0)
//            {
//              elem->getElementAssembledQty(FESystemElem::STRUCTURAL_F_WEIGHT_VECTOR::num(), 
//                                           &vector);
//              vec_it->second += vector;
//            }
          
        }
          break;
          
        case FORCE_VECTOR_SENSITIVITY_ENUM_ID:
        {
          if (nodal_loads.count(NODAL_TEMPERATURE::num()) > 0)
            {
              elem->getElementAssembledQtySensitivity
              (FESystemElem::STRUCTURAL_F_T_VECTOR::num(), &vector);
              vec_it->second += vector;
            }
          
          if (surface_loads.count(SURFACE_PRESSURE::num()) > 0)
            {
              elem->getElementAssembledQtySensitivity
              (FESystemElem::STRUCTURAL_F_PRESSURE_VECTOR::num(), &vector);
              vec_it->second += vector;
            }

//          if (volume_loads.count(SELF_WEIGHT::num()) > 0)
//            {
//              elem->getElementAssembledQtySensitivity
//              (FESystemElem::STRUCTURAL_F_WEIGHT_VECTOR::num(), &vector);
//              vec_it->second += vector;
//            }
        }
          break;
          
        default:
          abort();
          break;
      }
    }
  
  static std::map<unsigned int, double >::iterator real_it, real_end;
  real_it = real_elem_data.begin();
  real_end = real_elem_data.end();
  
  for (; real_it != real_end; real_it++)
    {
      real_it->second = 0.0;
      
      switch(real_it->first)
      {
        case MODEL_MASS_ENUM_ID:
          real_it->second += elem->getElementMass(false);
          break;
          
        case MODEL_MASS_SENSITIVITY_ENUM_ID:
          real_it->second += elem->getElementMass(true);
          break;
          
        default:
          abort();
      }
    }      
  
  elem->clearLoadAndDofs();
}




void
Discipline::StructuralAnalysis::performAdditionalOperationsOnGlobalQuantities
(std::map<unsigned int, SparseMatrix<double>* >& matrix_map,
 std::map<unsigned int, NumericVector<double>* >& vector_map)
{
  std::map<unsigned int, SparseMatrix<double>* >::iterator matrix_it, matrix_end;  
  matrix_it = matrix_map.begin();
  matrix_end = matrix_map.end();
  
  
  // iterate over the quantities and perform the operations on them
  for ( ; matrix_it != matrix_end; matrix_it++)
    {
      switch (matrix_it->first)
      {
        case TRANSIENT_C2_MATRIX_ENUM_ID:
        case SYSTEM_MATRIX_ENUM_ID:
        case JACOBIAN_MATRIX_ENUM_ID:
        case SYSTEM_MATRIX_SENSITIVITY_ENUM_ID:
        case EIGENPROBLEM_A_MATRIX_ENUM_ID:
        case EIGENPROBLEM_B_MATRIX_ENUM_ID:
        case EIGENPROBLEM_A_MATRIX_SENSITIVITY_ENUM_ID:
        case EIGENPROBLEM_B_MATRIX_SENSITIVITY_ENUM_ID:
        {
          // keep going
        }
          break;
          
        default:
          Assert(false, 
                 FESystemExceptions::ExcEnumCaseNotImplemented
                 (Driver::AnalysisDriverQtyEnum::enumName(matrix_it->first)));
      }
    }
  
  
  
  std::map<unsigned int, NumericVector<double>* >::iterator vector_it, vector_end;
  vector_it = vector_map.begin();
  vector_end = vector_map.end();
  
  for ( ; vector_it != vector_end; vector_it++)
    {
      switch (vector_it->first)
      {
        case FORCE_VECTOR_ENUM_ID:
        {
          Loads::NodalLoadDataInfo nodal_load_info;
          
          nodal_load_info.setLoadCaseID(this->analysis_driver->getCurrentLoadCase());
          nodal_load_info.setLoadNameEnumID(NODAL_FORCE::num());
          nodal_load_info.setLoadClassEnumID(NODAL_POINT_LOAD::num());
          nodal_load_info.setNDofs(6);
          if (this->solution_base->getDisciplineTransientNatureEnumID
              (Discipline::STRUCTURAL_DISCIPLINE::num()) == Discipline::TRANSIENT_SOLUTION::num())
            {
              Driver::TransientAnalysisDriver* transient_driver = 
              dynamic_cast<Driver::TransientAnalysisDriver*>(this->analysis_driver);
              nodal_load_info.setTime(transient_driver->getCurrentAnalysisTime());
            }
          
          FESystemUtility::AutoPtrVector<Loads::NodalLoadCombination> 
          point_loads(this->fesystem_controller.load_database->getAllLoadCombinations
                      <Loads::NodalLoadDataInfo, Loads::NodalLoadCombination>
                      (nodal_load_info).release());
          
          unsigned int node =0, dof =0; double value = 0.0;
          
          std::vector<Loads::NodalLoadCombination*>::const_iterator load_it, load_end;
          load_it = point_loads.get()->begin();
          load_end = point_loads.get()->end();
          
          for (; load_it != load_end ; load_it++)
            {
              node =  (**load_it).getNodeID();
              
              // now, iterate for each dof
              for (unsigned int dof_it=0; dof_it < 6; dof_it++)
                {
                  // get the load values and the dof index in the global load vector
                  // for this load
                  value = dynamic_cast<Loads::NodalPointLoadCombination*>(*load_it)->getValue(dof_it);
                  dof = (this->mesh_data->getNodeFromForeignID(node))->dof_number(0,dof_it,0);
                  
                  // add the calculated load to the global load_vector 
                  vector_it->second->add(dof, value);	
                }
            }
        }
          break;
          
          
        case FORCE_VECTOR_SENSITIVITY_ENUM_ID:
        {
          Loads::NodalLoadDataInfo nodal_load_info;
          
          nodal_load_info.setLoadCaseID(this->analysis_driver->getCurrentLoadCase());
          nodal_load_info.setLoadNameEnumID(NODAL_FORCE::num());
          nodal_load_info.setLoadClassEnumID(NODAL_POINT_LOAD::num());
          nodal_load_info.setNDofs(6);
          nodal_load_info.setDVID(this->analysis_driver->getCurrentDesignParameter().getID());
          if (this->solution_base->getDisciplineTransientNatureEnumID
              (Discipline::STRUCTURAL_DISCIPLINE::num()) == Discipline::TRANSIENT_SOLUTION::num())
            {
              Driver::TransientAnalysisDriver* transient_driver = 
              dynamic_cast<Driver::TransientAnalysisDriver*>(this->analysis_driver);
              nodal_load_info.setTime(transient_driver->getCurrentAnalysisTime());
            }
          
          FESystemUtility::AutoPtrVector<Loads::NodalLoadCombination> 
          point_loads(this->fesystem_controller.load_database->getAllLoadCombinations
                      <Loads::NodalLoadDataInfo, Loads::NodalLoadCombination>
                      (nodal_load_info).release());
          
          unsigned int node =0, dof =0; double value = 0.0;
          
          std::vector<Loads::NodalLoadCombination*>::const_iterator load_it, load_end;
          load_it = point_loads.get()->begin();
          load_end = point_loads.get()->end();
          
          for (; load_it != load_end ; load_it++)
            {
              node =  (**load_it).getNodeID();
              
              // now, iterate for each dof
              for (unsigned int dof_it=0; dof_it < 6; dof_it++)
                {
                  // get the load values and the dof index in the global load vector
                  // for this load
                  value = dynamic_cast<Loads::NodalPointLoadCombination*>(*load_it)->getValue(dof_it);
                  dof = (this->mesh_data->getNodeFromForeignID(node))->dof_number(0,dof_it,0);
                  
                  // add the calculated load to the global load_vector 
                  vector_it->second->add(dof, value);
                }
            }
        }
          break;
          
        default:
          Assert(false, 
                 FESystemExceptions::ExcEnumCaseNotImplemented
                 (Driver::AnalysisDriverQtyEnum::enumName(matrix_it->first)));
      }
    }
}



const NumericVector<double>& 
Discipline::StructuralAnalysis::getSolutionVector(const unsigned int transient_order,
                                                  const bool sensitivity)
{
  if (transient_order > 0)
    Assert(this->solution_base->getDisciplineTransientNatureEnumID
           (Discipline::STRUCTURAL_DISCIPLINE::num()) == 
           Discipline::TRANSIENT_SOLUTION::num(), 
           ExcInternalError());
  
  // a pointer to the data info object
  static FESystemDatabase::DataInfoBase *data_info_ptr = NULL;
  
  // check the current operation being performed
  switch (this->analysis_driver->getCurrentAnalysisKind())
  {
    case ANALYSIS_ENUM_ID:
    {
      // sensitivity information for analysis should not be needed for any solution
      Assert(!sensitivity, ExcInternalError());
      
      switch (this->solution_base->getDisciplineTransientNatureEnumID
              (Discipline::STRUCTURAL_DISCIPLINE::num()))
      {
        case STEADY_STATE_SOLUTION_ENUM_ID:
        {
          static FESystemDatabase::TimeIndependentDataInfo data_info;
          data_info.clear();
          data_info_ptr = &data_info;
          switch (this->solution_base->getSolutionEnumID())
          {
            case STRUCTURAL_VIBRATION_EIGEN_SOLUTION_ENUM_ID:
            case LINEARIZED_BUCKLING_EIGEN_SOLUTION_ENUM_ID:
            case AEROELASTICITY_SOLUTION_ENUM_ID:
            {
              // this solution should have been performed, and hence, its
              // name should be available from the solution
              data_info.setName("Solution");
            }
              break;
              
            default:
              Assert(false, ExcInternalError());
          }
        }
          break;
          
        case TRANSIENT_SOLUTION_ENUM_ID:
        default:
          Assert(false, ExcInternalError());
      }
    }
      break;
      
    case SENSITIVITY_ANALYSIS_ENUM_ID:
    {
      switch (this->solution_base->getDisciplineTransientNatureEnumID
              (Discipline::STRUCTURAL_DISCIPLINE::num()))
      {
        case STEADY_STATE_SOLUTION_ENUM_ID:
        {
          static FESystemDatabase::TimeIndependentDataInfo data_info;
          data_info.clear();
          data_info_ptr = &data_info;
          switch (this->solution_base->getSolutionEnumID())
          {
            case LINEAR_STRESS_SOLUTION_ENUM_ID:
            case STRUCTURAL_VIBRATION_EIGEN_SOLUTION_ENUM_ID:
            case LINEARIZED_BUCKLING_EIGEN_SOLUTION_ENUM_ID:
            case AEROELASTICITY_SOLUTION_ENUM_ID:
            {
              data_info.setName("Solution");
            }
              break;
              
            default:
              Assert(false, ExcInternalError());
          }
        }
          break;
          
        case TRANSIENT_SOLUTION_ENUM_ID:
        default:
          Assert(false, ExcInternalError());
      }
    }
      break;
      
    case POST_PROCESS_ENUM_ID:
    {
      // needs to be implemented
      Assert(false, ExcInternalError());
    }
      break;
      
    default:
      Assert(false, ExcInternalError());
  }
  
  if (sensitivity)
    data_info_ptr->setDVID(this->analysis_driver->getCurrentDesignParameter().getID());
  
  data_info_ptr->setLoadCase(this->analysis_driver->getCurrentLoadCase());
  data_info_ptr->setDisciplineEnumID(Discipline::STRUCTURAL_DISCIPLINE::num());
  return this->loaded_solution.loadSolutionFromDatabase
  (*(this->fesystem_controller.global_data_storage.get()), *data_info_ptr);
}



bool
Discipline::StructuralAnalysis::disciplineHasMatrix(const unsigned int qty_enum_ID) const
{
  bool return_val = false;
  
  switch (qty_enum_ID)
  {
    case TRANSIENT_C2_MATRIX_ENUM_ID:
    {
      switch (this->solution_base->getSolutionEnumID())
      {
        case TRANSIENT_STRUCTURAL_SOLUTION_ENUM_ID:
        case AEROELASTICITY_SOLUTION_ENUM_ID:
          return_val = true;
          break;
          
        default:
          break;
      }
    }
      break;
      
    case TRANSIENT_C1_MATRIX_ENUM_ID:
      // for now, no damping is provided in the structure
      return_val = false;
      break;
      
    case SYSTEM_MATRIX_ENUM_ID:
      return_val = true;
      break;
      
    default:
      Assert(false, ExcInternalError());
  }
  
  return return_val;
}






std::auto_ptr<FESystemDatabase::DataInfoBase> 
Discipline::StructuralAnalysis::getDataInfoForQty(const unsigned int qty_enum_ID)
{
  std::auto_ptr<FESystemDatabase::DataInfoBase> data_info;
  bool set_transient_data = false,
  set_DV_ID = false;
  
  switch (qty_enum_ID)
  {
    case TRANSIENT_C2_MATRIX_ENUM_ID:
    {
      // check the kind of analysis in progress
      switch (this->solution_base->getDisciplineTransientNatureEnumID
              (Discipline::STRUCTURAL_DISCIPLINE::num()))
      {
        case STEADY_STATE_SOLUTION_ENUM_ID:
        {
          data_info.reset(new FESystemDatabase::TimeIndependentDataInfo);
          data_info->setName("M_Matrix");
        }
          break;
          
        case TRANSIENT_SOLUTION_ENUM_ID:
        {
          FESystemDatabase::TimeDependentDataInfo *ptr = new FESystemDatabase::TimeDependentDataInfo;
          data_info.reset(ptr);
          data_info->setName("M_Matrix");
          ptr->setOrder(0);
          set_transient_data = true;
          ptr = NULL;
        }
          break;
          
        default:
          // at this point of time, no other analysis needs the specification of the elem
          // dof vector
          Assert(false, ExcInternalError());
          
      }
    }
      break;
      
    case TRANSIENT_C2_MATRIX_SENSITIVITY_ENUM_ID:
    {
      data_info.reset(this->getDataInfoForQty(Driver::TRANSIENT_C2_MATRIX::num()).release());
      set_DV_ID = true;
    }
      break;
      
    case SYSTEM_MATRIX_ENUM_ID:
    {
      // check the kind of analysis in progress
      switch (this->solution_base->getDisciplineTransientNatureEnumID
              (Discipline::STRUCTURAL_DISCIPLINE::num()))
      {
        case STEADY_STATE_SOLUTION_ENUM_ID:
        {
          data_info.reset(new FESystemDatabase::TimeIndependentDataInfo);
          data_info->setName("K_Matrix");
        }
          break;
          
        case TRANSIENT_SOLUTION_ENUM_ID:
        {
          FESystemDatabase::TimeDependentDataInfo *ptr = new FESystemDatabase::TimeDependentDataInfo;
          data_info.reset(ptr);
          data_info->setName("K_Matrix");
          ptr->setOrder(0);
          set_transient_data = true;
          ptr = NULL;
        }
          break;
          
        default:
          // at this point of time, no other analysis needs the specification of the elem
          // dof vector
          Assert(false, ExcInternalError());
          
      }
    }
      break;
      
      
    case SYSTEM_MATRIX_SENSITIVITY_ENUM_ID:
    {
      data_info.reset(this->getDataInfoForQty(Driver::SYSTEM_MATRIX::num()).release());
      set_DV_ID = true;
    }
      break;
      
      
      
    case EIGENPROBLEM_A_MATRIX_ENUM_ID:
    {
      switch (this->solution_base->getSolutionEnumID())
      {
        case STRUCTURAL_VIBRATION_EIGEN_SOLUTION_ENUM_ID:
        case LINEARIZED_BUCKLING_EIGEN_SOLUTION_ENUM_ID:
        {
          data_info.reset(new FESystemDatabase::TimeIndependentDataInfo);
          data_info->setName("K_Matrix");
        }
          break;
          
        default:
          Assert(false, ExcInternalError());
      }
    }
      break;
      
      
      
    case EIGENPROBLEM_A_MATRIX_SENSITIVITY_ENUM_ID:
    {
      data_info.reset(this->getDataInfoForQty(Driver::EIGENPROBLEM_A_MATRIX::num()).release());
      set_DV_ID = true;
    }
      break;
      
      
    case EIGENPROBLEM_B_MATRIX_ENUM_ID:
    {
      switch (this->solution_base->getSolutionEnumID())
      {
        case STRUCTURAL_VIBRATION_EIGEN_SOLUTION_ENUM_ID:
        {
          data_info.reset(new FESystemDatabase::TimeIndependentDataInfo);
          data_info->setName("Mass_Matrix");
        }
          break;
          
        case LINEARIZED_BUCKLING_EIGEN_SOLUTION_ENUM_ID:
        {
          data_info.reset(new FESystemDatabase::TimeIndependentDataInfo);
          data_info->setName("K_G_Matrix");
        }
          break;
          
        default:
          Assert(false, ExcInternalError());
      }
    }
      break;
      
      
      
    case EIGENPROBLEM_B_MATRIX_SENSITIVITY_ENUM_ID:
    {
      data_info.reset(this->getDataInfoForQty(Driver::EIGENPROBLEM_B_MATRIX::num()).release());
      set_DV_ID = true;
    }
      break;
      
      
      
    case FORCE_VECTOR_ENUM_ID:
    {
      // check the kind of analysis in progress
      switch (this->solution_base->getDisciplineTransientNatureEnumID
              (Discipline::STRUCTURAL_DISCIPLINE::num()))
      {
        case STEADY_STATE_SOLUTION_ENUM_ID:
        {
          data_info.reset(new FESystemDatabase::TimeIndependentDataInfo);
          data_info->setName("F_Vector");
        }
          break;
          
        case TRANSIENT_SOLUTION_ENUM_ID:
        {
          FESystemDatabase::TimeDependentDataInfo *ptr = new FESystemDatabase::TimeDependentDataInfo;
          data_info.reset(ptr);
          data_info->setName("F_Vector");
          ptr->setOrder(0);
          set_transient_data = true;
          ptr = NULL;
        }
          break;
          
        default:
          // at this point of time, no other analysis needs the specification of the elem
          // dof vector
          Assert(false, ExcInternalError());
          
      }
    }
      break;
      
    case FORCE_VECTOR_SENSITIVITY_ENUM_ID:
    {
      data_info.reset(this->getDataInfoForQty(Driver::FORCE_VECTOR::num()).release());
      set_DV_ID = true;
    }
      break;
			
			
    default:
      Assert(false, ExcInternalError());
      break;
  }
  
  
  if (set_DV_ID)
    data_info->setDVID(this->analysis_driver->getCurrentDesignParameter().getID());
  else
    {
      data_info->setDisciplineEnumID(Discipline::STRUCTURAL_DISCIPLINE::num());
      data_info->setLoadCase(this->analysis_driver->getCurrentLoadCase());
    }
  
  if (set_transient_data)
    {
      Driver::TransientAnalysisDriver* transient_driver = 
      dynamic_cast<Driver::TransientAnalysisDriver*>(this->analysis_driver);
      FESystemDatabase::TimeDependentDataInfo* time_dep_info = 
      dynamic_cast<FESystemDatabase::TimeDependentDataInfo*>(data_info.get());
      time_dep_info->setTransientIterationInfo(transient_driver->getCurrentIterationNumber(),
                                               transient_driver->getCurrentAnalysisTime());
    }
  
  
  return data_info;
}



void
Discipline::StructuralAnalysis::
setMeshNodalIncrementsFromDisplacementSolution(const NumericVector<double>& vec)
{
  // iterate over each node, get the u,v,w values from the given solution vector and 
  // use them to update the x,y,z locations of the nodes
  MeshBase::node_iterator it = this->mesh->nodes_begin();
  MeshBase::node_iterator end = this->mesh->nodes_end();

  unsigned int dof_num = 0;
  
  for ( ; it != end; it++)
    {
      for (unsigned int dof_it=0; dof_it < 3; dof_it++)
        {
          dof_num = (**it).dof_number(0,dof_it,0);
          (**it)(dof_it) += vec(dof_num);
        }
    }
}



unsigned int
Discipline::StructuralAnalysis::getTransientSystemOrder()
{
  unsigned int order = 0;
  // return the value based on the solution being performed
  switch (this->solution_base->getSolutionEnumID()) 
  {
    case TRANSIENT_SOLUTION_ENUM_ID:
      order = 2;
      break;
      
    case STRUCTURAL_VIBRATION_EIGEN_SOLUTION_ENUM_ID:
    case LINEARIZED_BUCKLING_EIGEN_SOLUTION_ENUM_ID:
    case STEADY_STATE_SOLUTION_ENUM_ID:
    case NONLINEAR_STATIC_STRUCTURAL_SOLUTION_ENUM_ID:
    case AEROELASTICITY_SOLUTION_ENUM_ID:
      order = 0;
      break;
      
    default:
      Assert (false, ExcInternalError());
      break;
  }
  return order;
}




bool 
Discipline::StructuralAnalysis::sameSystemMatrixForLoadCases()
{
  switch (this->analysis_type_enum_ID)
  {
    case LINEAR_ANALYSIS_ENUM_ID:
    {
      if (this->checkPropertyDependenceOnTemperature())
        return false;
      else
        return true;
    }
      break;
      
    case NONLINEAR_ANALYSIS_ENUM_ID:
      return false;
      break;
      
    default:
      abort();
      break;
  }
}

