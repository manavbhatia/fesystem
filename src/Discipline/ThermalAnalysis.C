// $Id: ThermalAnalysis.C,v 1.17.4.9 2008-06-03 05:19:33 manav Exp $

// C++ includes
#include <sstream>

// FESystem includes
#include "Discipline/ThermalAnalysis.h"
#include "Discipline/ThermalDisciplineInfo.h"
#include "FESystem/FESystemController.h"
#include "FESystem/AnalysisCase.h"
#include "Mesh/MeshList.h"
#include "Mesh/FEMeshData.h"
#include "Database/GlobalDataStorage.h"
#include "AnalysisDriver/AnalysisDriver.h"
#include "AnalysisDriver/LinearAnalysisDriver.h"
#include "AnalysisDriver/NonLinearAnalysisDriver.h"
#include "AnalysisDriver/TransientAnalysisDriver.h"
#include "Radiation/RadiationCavityAnalysis.h"
#include "Radiation/RadiationCavity.h"
#include "ThermalElems/thermal_elem.h"
#include "FESystem/FESystemExceptions.h"
#include "Loads/LoadDatabase.h"
#include "Loads/LoadDataInfo.h"
#include "Solutions/SolutionBase.h"
#include "Numerics/PetscSeqVector.h"
#include "Numerics/PetscSeqDenseMatrix.h"

#include "Loads/LoadCombination.h"

// libMesh includes
#include "geom/edge_edge2.h"
#include "numerics/sparse_matrix.h"
#include "numerics/numeric_vector.h"
#include "numerics/dense_matrix.h"
#include "numerics/dense_vector.h"
#include "numerics/petsc_matrix.h"


Discipline::ThermalAnalysis::ThermalAnalysis(FESystem::FESystemController& controller,
                                             const Discipline::ThermalDisciplineInfo& info):
Discipline::AnalysisDisciplineBase(controller, info)
{
  
  // if this analysis is nonlinear, init the vector of radiation cavity analysis
  if (this->analysis_type_enum_ID == Discipline::NONLINEAR_ANALYSIS::num())
    {
      
      // get the vector of radiation cavities in the analysis 
      const std::vector<unsigned int>& radiation_cavity_IDs = info.getRadiationCavityIDs();
      
      std::vector<unsigned int>::const_iterator it, end;
      it = radiation_cavity_IDs.begin();
      end = radiation_cavity_IDs.end();
      
      for (; it != end; it++)
        {
          RadiationCavity& cavity = 
          *(this->fesystem_controller.mesh_list->getRadiationCavity(*it));
          
          std::auto_ptr<RadiationCavityAnalysis> 
          cavity_analysis(new RadiationCavityAnalysis(this->fesystem_controller, *this, cavity));
          
          this->radiation_cavity_analysis_vector.push_back(cavity_analysis.release());
          
          // get the node set list for each cavity, and add rod elements connecting 
          // each one of those nodes
          // this is done to take care of the connectivity between nodes 
          // to take care of the radiation exchange between them
          
          const std::vector<Node*>& radiation_nodes = cavity.getFENodes();
          
          std::vector<Node*>::const_iterator node_it, node_end, nested_node_it;
          node_it = radiation_nodes.begin();
          node_end = radiation_nodes.end();
          
          Elem* this_elem = NULL;
          bool insert_return = false;
          
          for (; node_it != node_end; node_it++)
            {
              nested_node_it = node_it;
              nested_node_it++;
              for ( ; nested_node_it != node_end; nested_node_it++)
                {
                  // create an elem connecting these two nodes
                  this_elem = this->mesh->add_elem(new Edge2);
                  this_elem->set_node(0) = *node_it;
                  this_elem->set_node(1) = *nested_node_it;
                  
                  // add this elem to the dummy element set
                  insert_return = this->dummy_element_set->insert(this_elem).second;
                  assert (insert_return == true);
                }
            }
        }
    }
  
  this->addAnalysisVariablesAndInitialize();
}





Discipline::ThermalAnalysis::~ThermalAnalysis()
{
  // iterate over the radiation cavity analysis and delete them
  std::vector<RadiationCavityAnalysis*>::iterator it, end;
  it = this->radiation_cavity_analysis_vector.begin();
  end = this->radiation_cavity_analysis_vector.end();
  
  for (; it != end; it++)
    {
      delete *it;
      *it = NULL;
    }
}





void 
Discipline::ThermalAnalysis::addAnalysisVariables()
{  
  FEType fetype(FIRST, LAGRANGE);
  
  // add a temperature variable to the analysis driver
  this->addVariable("Temperature", fetype);
}






unsigned int
Discipline::ThermalAnalysis::getEigenProblemKindEnumID()
{
  // this is not pertinent for thermal discipline, since no eigenproblem has yet 
  // been defined for this. Hence, it is an error to call this method for thermal analysis
  Assert(false, ExcInternalError());
  return 0;
}




std::auto_ptr<Loads::DirichletBoundaryConditionDataInfo> 
Discipline::ThermalAnalysis::getBoundaryConditionLoadInfo() const
{
  Assert(this->analysis_driver != NULL, ExcInternalError());
  
  std::auto_ptr<Loads::DirichletBoundaryConditionDataInfo> return_ptr;
  
  return_ptr.reset(new Loads::DirichletBoundaryConditionDataInfo);
  return_ptr->setLoadCaseID(this->analysis_driver->getCurrentLoadCase());
  return_ptr->setLoadNameEnumID(TEMPERATURE_BOUNDARY_CONDITION::num());
  return_ptr->setLoadClassEnumID(DIRICHLET_BOUNDARY_CONDITION::num());
  
  // check if the analysis is time dependent
  switch (this->solution_base->getDisciplineTransientNatureEnumID
          (Discipline::THERMAL_DISCIPLINE::num()))
  {
    case STEADY_STATE_SOLUTION_ENUM_ID:
    {
      //nothing to be done here
    }
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
Discipline::ThermalAnalysis::calculateTransientMatrix(const unsigned int order, bool sensitivity)
{
  // this is to avoid getting warnings about unused parameters
  (void) sensitivity;
  (void) order;
  // for now, everything will be calculated. the below logic does not take into account 
  // nonlinear analysis. that will be done later.
  return true;
  
  // needs to be implemented to allow for some fine grained control over when to 
  // reclaculated the quantity
  //   unsigned int load_case = this->analysis_driver->getCurrentLoadCase();
	
  //   // get the name by which the quantity will be referred
  //   std::string name = this->getStringNameForQty(Driver::SYSTEM_MATRIX::num(), sensitivity);
	
  //   // this bool will check if this quantity needs to be recalculated and resaved
  //   // if the quantity is dependent on temperature, and a sensitivity analysis 
  //   // is being performed, then it is assumed that the quantity already exists, and was 
  //   // calculated at the converged temperature for that load case. 
  //   bool recalculate_and_store = false;
  //   if (this->checkPropertyDependenceOnTemperature())
  //     if (!sensitivity && 
  //         this->analysis_driver->getCurrentAnalysisKind() != Driver::SENSITIVITY_ANALYSIS::num())
  //       recalculate_and_store = true;
  
  
  //   // check if the K_C + K_h matrix exists in the database
  //   bool K_c_K_h_matrix_present = 
  //     this->fesystem_controller.global_data_storage->checkIfMatrixExists
  //     (load_case, name);
  
  //   // if not, calculate and save it
  //   if (! K_c_K_h_matrix_present || recalculate_and_store)
  //     return true;
  //   else 
  //     return false;
}




bool 
Discipline::ThermalAnalysis::calculateK(bool sensitivity)
{
  // this is to avoid getting warnings about unused parameters
  (void) sensitivity;
  // for now, everything will be calculated. the below logic does not take into account 
  // nonlinear analysis. that will be done later.
  return true;
  
  // needs to be implemented to allow for some fine grained control over when to 
  // reclaculated the quantity
  //   unsigned int load_case = this->analysis_driver->getCurrentLoadCase();
	
  //   // get the name by which the quantity will be referred
  //   std::string name = this->getStringNameForQty(Driver::SYSTEM_MATRIX::num(), sensitivity);
	
  //   // this bool will check if this quantity needs to be recalculated and resaved
  //   // if the quantity is dependent on temperature, and a sensitivity analysis 
  //   // is being performed, then it is assumed that the quantity already exists, and was 
  //   // calculated at the converged temperature for that load case. 
  //   bool recalculate_and_store = false;
  //   if (this->checkPropertyDependenceOnTemperature())
  //     if (!sensitivity && 
  //         this->analysis_driver->getCurrentAnalysisKind() != Driver::SENSITIVITY_ANALYSIS::num())
  //       recalculate_and_store = true;
  
  
  //   // check if the K_C + K_h matrix exists in the database
  //   bool K_c_K_h_matrix_present = 
  //     this->fesystem_controller.global_data_storage->checkIfMatrixExists
  //     (load_case, name);
  
  //   // if not, calculate and save it
  //   if (! K_c_K_h_matrix_present || recalculate_and_store)
  //     return true;
  //   else 
  //     return false;
}





bool 
Discipline::ThermalAnalysis::calculateJac()
{
  // for now, everything will be calculated. the below logic does not take into account 
  // nonlinear analysis. that will be done later.
  return true;
  
  //   // get the name by which the quantity will be referred
  //   std::string name = this->getStringNameForQty(Driver::JACOBIAN_MATRIX::num(), false);
	
  //   // check if the jacobian matrix exists in the database for this iteration
  //   // this bool will check if this quantity needs to be recalculated and resaved
  //   // if the quantity is dependent on temperature, and a sensitivity analysis 
  //   // is being performed, then it is assumed that the quantity already exists, and was 
  //   // calculated at the converged temperature for that load case. 
  //   bool recalculate_and_store = false;
  //   if (this->analysis_driver->getCurrentAnalysisKind() != Driver::SENSITIVITY_ANALYSIS::num())
  //     recalculate_and_store = true;
  
  //   return recalculate_and_store;
}




bool 
Discipline::ThermalAnalysis::calculateEigenProblemAMatrix(bool sensitivity)
{
  // this is to avoid warnings about unused parameter.
  (void) sensitivity;
  // there still is no nonlinear structural analysis capability yet.
  abort();
}



bool
Discipline::ThermalAnalysis::calculateEigenProblemBMatrix(bool sensitivity)
{
  // this is to avoid warnings about unused parameter.
  (void) sensitivity;
  // there still is no nonlinear structural analysis capability yet.
  abort();
}



bool 
Discipline::ThermalAnalysis::calculateF(bool sensitivity)
{
  // this is to avoid warnings about unused parameter.
  (void) sensitivity;  
  // for now, everything will be calculated. the below logic does not take into account 
  // nonlinear analysis. that will be done later.
  return true;
  
  // this needs to be reimplemented to allow for some fine grained control over
  // when to calculate the quantity
  //   unsigned int load_case = this->analysis_driver->getCurrentLoadCase();
	
  //   // get the name by which the quantity will be referred
  //   std::string name = this->getStringNameForQty(Driver::FORCE_VECTOR::num(), sensitivity);
	
  //   // this bool will check if this quantity needs to be recalculated and resaved
  //   // if the quantity is dependent on temperature, and a sensitivity analysis 
  //   // is being performed, then it is assumed that the quantity already exists, and was 
  //   // calculated at the converged temperature for that load case. 
  //   bool recalculate_and_store = false;
  //   if (this->checkPropertyDependenceOnTemperature())
  //     if (!sensitivity && 
  //         this->analysis_driver->getCurrentAnalysisKind() != Driver::SENSITIVITY_ANALYSIS::num())
  //       recalculate_and_store = true;
  
  //   // check if the vector exists in the database
  //   bool Force_vector_present = 
  //     this->fesystem_controller.global_data_storage->checkIfVectorExists
  //     (load_case, name);
  
  //   // if not, calculate and save it
  //   if (! Force_vector_present || recalculate_and_store)
  //     return true;
  //   else 
  //     return false;
}




void 
Discipline::ThermalAnalysis::getElemQty
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
  static DenseVector<double> vector, elem_dofs, elem_dof_sensitivity, elem_dofs_deriv;
  static std::vector<DenseVector<double>*> elem_dof_vec(2), elem_dof_sens_vec(2);
  
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
      volume_load_map[VOLUME_HEAT_LOAD::num()] = SCALAR_VOLUME_LOAD::num();
      surface_load_map[SURFACE_HEAT_LOAD::num()] = SCALAR_SURFACE_LOAD::num();
      surface_load_map[SURFACE_RADIATION_HEAT_LOAD::num()] = SURFACE_RADIATION_LOAD::num();
      surface_load_map[SURFACE_CONVECTION_HEAT_LOAD::num()] = SURFACE_CONVECTION_LOAD::num();
      
      created_load_maps = true;
    }
  
  elem_ID = elem->getID();
  load_case = this->analysis_driver->getCurrentLoadCase();
  
  // set the transient nature of this discipline
  static double time_val = 0.0;
  static bool transient_nature = false;
  
  if (this->solution_base->getDisciplineTransientNatureEnumID(Discipline::THERMAL_DISCIPLINE::num()) 
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
                            &volume_load_sensitivity, &surface_load_sensitivity, &nodal_load_sensitivity);
    }
  
  if (this->analysis_type_enum_ID == NONLINEAR_ANALYSIS::num())
    {
      // all quantities will be evaluated at the same dof value.
      this->getElemDofValues(elem->getElem(FESystemElem::BASE_ELEM::num()),
                             elem_dofs, 0, false);
      
      elem_dof_vec[0] = &elem_dofs;
      
      if (this->solution_base->getDisciplineTransientNatureEnumID
          (Discipline::THERMAL_DISCIPLINE::num()) == 
          Discipline::TRANSIENT_SOLUTION::num())
        {
          // all quantities will be evaluated at the same dof value.
          this->getElemDofValues(elem->getElem(FESystemElem::BASE_ELEM::num()),
                                 elem_dofs_deriv, 1, false);
          elem_dof_vec[1] = &elem_dofs_deriv;
        }
      else
        elem_dof_vec[1] = NULL;
      
      elem_dof_sens_vec[0] = NULL;
      elem_dof_sens_vec[1] = NULL;
      
      elem->setDofValues(elem_dof_vec, elem_dof_sens_vec);
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
        case TRANSIENT_C1_MATRIX_ENUM_ID:
        {
          // get K_c and K_h
          elem->getElementAssembledQty(FESystemElem::THERMAL_C_MATRIX::num(), &matrix);
          mat_it->second += matrix;
        }
          break;
          
        case TRANSIENT_C1_MATRIX_SENSITIVITY_ENUM_ID:
        {
          // get K_c and K_h
          elem->getElementAssembledQtySensitivity(FESystemElem::THERMAL_C_MATRIX::num(), &matrix);
          mat_it->second += matrix;
        }
          break;
          
          
        case SYSTEM_MATRIX_ENUM_ID:
        {
          // get K_c and K_h
          elem->getElementAssembledQty(FESystemElem::THERMAL_K_C_MATRIX::num(), &matrix);
          mat_it->second += matrix;
          
          if (surface_loads.count(SURFACE_CONVECTION_HEAT_LOAD::num()) > 0)
            {
              elem->getElementAssembledQty(FESystemElem::THERMAL_K_H_MATRIX::num(), &matrix);
              mat_it->second += matrix;
            }
        }
          break;
          
          
        case JACOBIAN_MATRIX_ENUM_ID:
        {
          // get K_c and K_h
          elem->getElementAssembledQty(FESystemElem::THERMAL_K_C_MATRIX::num(), &matrix);
          mat_it->second += matrix;
          
          if (surface_loads.count(SURFACE_CONVECTION_HEAT_LOAD::num()) > 0)
            {
              elem->getElementAssembledQty(FESystemElem::THERMAL_K_H_MATRIX::num(), &matrix);
              mat_it->second += matrix;
            }
          
          // get F_sigma_Jac 
          if (surface_loads.count(SURFACE_RADIATION_HEAT_LOAD::num()) > 0)
            {
              elem->getElementAssembledQty(FESystemElem::THERMAL_EMITTED_RAD_JAC_MATRIX::num(),
                                           &matrix);
              matrix.scale(-1.0);
              mat_it->second += matrix;
            }
          
          if (this->checkPropertyDependenceOnTemperature())
            {
              // get the K_c jacobian due to temperature dependence of material property
              elem->getElementAssembledQty(FESystemElem::THERMAL_K_C_JAC_MATRIX::num(), &matrix);
              mat_it->second += matrix;
              
              if (this->solution_base->getDisciplineTransientNatureEnumID
                  (Discipline::THERMAL_DISCIPLINE::num()) == 
                  Discipline::TRANSIENT_SOLUTION::num() &&
                  this->analysis_type_enum_ID == NONLINEAR_ANALYSIS::num())
                {
                  // get the K_c jacobian due to temperature dependence of material property
                  elem->getElementAssembledQty(FESystemElem::THERMAL_C_JAC_MATRIX::num(), &matrix);
                  mat_it->second += matrix;
                }
              
            }
        }
          break;
          
        case SYSTEM_MATRIX_SENSITIVITY_ENUM_ID:
        {
          // get K_c and K_h sensitivity
          elem->getElementAssembledQtySensitivity(FESystemElem::THERMAL_K_C_MATRIX::num(),
                                                  &matrix);
          mat_it->second += matrix;
          
          if (surface_loads.count(SURFACE_CONVECTION_HEAT_LOAD::num()) > 0)
            {
              elem->getElementAssembledQtySensitivity(FESystemElem::THERMAL_K_H_MATRIX::num(),
                                                      &matrix);
              mat_it->second += matrix;
            }
        }
          break;
          
          
        default:
          Assert(false, ExcInternalError());
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
          // get F_h, F_q, F_Q, F_nodal, F_sigma
          if (surface_loads.count(SURFACE_CONVECTION_HEAT_LOAD::num()) > 0)
            {
              elem->getElementAssembledQty(FESystemElem::THERMAL_F_H_VECTOR::num(), &vector);
              vec_it->second += vector;
            }
          
          if (volume_loads.count(VOLUME_HEAT_LOAD::num()) > 0)
            {
              elem->getElementAssembledQty(FESystemElem::THERMAL_F_VOL_VECTOR::num(), &vector);
              vec_it->second += vector;
            }
          
          if (surface_loads.count(SURFACE_HEAT_LOAD::num()) > 0)
            {
              elem->getElementAssembledQty(FESystemElem::THERMAL_F_Q_SURF_VECTOR::num(), &vector);
              vec_it->second += vector;
            }
          
          // this is only for nonlinear. Hence, an if block should be used to calculate
          // this only is nonlinear 
          // analysis is being performed
          if (this->analysis_type_enum_ID == NONLINEAR_ANALYSIS::num())
            if (surface_loads.count(SURFACE_RADIATION_HEAT_LOAD::num()) > 0)
              {
                elem->getElementAssembledQty(FESystemElem::THERMAL_F_EMITTED_RAD_VECTOR::num(),
                                             &vector);
                vec_it->second += vector;
              }
        }
          break;
          
        case FORCE_VECTOR_SENSITIVITY_ENUM_ID:
        {
          // get F_h, F_q, F_Q, F_nodal, F_sigma sensitivity
          if (surface_loads.count(SURFACE_CONVECTION_HEAT_LOAD::num()) > 0)
            {
              elem->getElementAssembledQtySensitivity(FESystemElem::THERMAL_F_H_VECTOR::num(),
                                                      &vector);
              vec_it->second += vector;
            }
          
          if (volume_loads.count(VOLUME_HEAT_LOAD::num()) > 0)
            {
              elem->getElementAssembledQtySensitivity(FESystemElem::THERMAL_F_VOL_VECTOR::num(),
                                                      &vector);
              vec_it->second += vector;
            }
          
          if (surface_loads.count(SURFACE_HEAT_LOAD::num()) > 0)
            {
              elem->getElementAssembledQtySensitivity(FESystemElem::THERMAL_F_Q_SURF_VECTOR::num(),
                                                      &vector);
              vec_it->second += vector;
            }
          
          if (this->analysis_type_enum_ID == NONLINEAR_ANALYSIS::num())
            if (surface_loads.count(SURFACE_RADIATION_HEAT_LOAD::num()) > 0)
              {
                elem->getElementAssembledQtySensitivity
                (FESystemElem::THERMAL_F_EMITTED_RAD_VECTOR::num(), &vector);
                vec_it->second += vector;
              }
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
Discipline::ThermalAnalysis::performAdditionalOperationsOnGlobalQuantities
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
        case JACOBIAN_MATRIX_ENUM_ID:
        {
          const std::vector<unsigned int>& nnz = this->dof_map->get_n_nz();
          std::vector<std::set<unsigned int> > ids(nnz.size());
          
          // get a reference to the solution vector
          NumericVector<double>& sol_vec = this->analysis_driver->getCurrentSolution(0);
          
          FESystemNumerics::PetscSeqVector<double> nodal_temp_vals;
          FESystemNumerics::PetscSeqDenseMatrix<double> radiation_jacobian;
          std::vector<unsigned int> rad_nodal_dofs;
          
          std::vector<RadiationCavityAnalysis*>::iterator it, end;
          it = this->radiation_cavity_analysis_vector.begin();
          end = this->radiation_cavity_analysis_vector.end();
          
          for (; it != end; it++)
            {
              // get the nodal dofs for the cavity
              rad_nodal_dofs.clear();
              (*it)->getFENodalDofIDs(rad_nodal_dofs);
              
              // prepare the nodal temperature vector for these dofs
              nodal_temp_vals.resize(rad_nodal_dofs.size());
              radiation_jacobian.resize(rad_nodal_dofs.size(),
                                        rad_nodal_dofs.size());
              
              for (unsigned int i=0; i < nodal_temp_vals.size(); i++)
                nodal_temp_vals.set(i, sol_vec(rad_nodal_dofs[i]));
              
              
              (*it)->getFEJacobianMatrix(nodal_temp_vals, radiation_jacobian);
              
              
              DenseMatrix<double> matrix(radiation_jacobian.m(),
                                         radiation_jacobian.n());
              
              for (unsigned int i=0; i < matrix.m(); i++)
                for (unsigned int j=0; j < matrix.n(); j++)
                  matrix(i,j) = radiation_jacobian.el(i,j);
              
              matrix_it->second->add_matrix(matrix, rad_nodal_dofs);
            }
        }
          break;
          
        default:
        {
          // keep going, nothing to be done otherwise
        }
          break;
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
          if (this->analysis_type_enum_ID == Discipline::NONLINEAR_ANALYSIS::num())
            {
              NumericVector<double>& sol_vec = this->analysis_driver->getCurrentSolution(0);
              FESystemNumerics::PetscSeqVector<double> radiation_load_vector, nodal_temp_vals;
              std::vector<unsigned int> rad_nodal_dofs;
              
              std::vector<RadiationCavityAnalysis*>::iterator it, end;
              it = this->radiation_cavity_analysis_vector.begin();
              end = this->radiation_cavity_analysis_vector.end();
              
              for (; it != end; it++)
                {
                  // get the nodal dofs for the cavity
                  rad_nodal_dofs.clear();
                  (*it)->getFENodalDofIDs(rad_nodal_dofs);
                  
                  // prepare the nodal temperature vector for these dofs
                  nodal_temp_vals.resize(rad_nodal_dofs.size());
                  radiation_load_vector.resize(rad_nodal_dofs.size());
                  
                  for (unsigned int i=0; i < nodal_temp_vals.size(); i++)
                    nodal_temp_vals.set(i, sol_vec(rad_nodal_dofs[i]));
                  
                  (*it)->getFENodalLoadVector(nodal_temp_vals,
                                              radiation_load_vector);
                  
                  // the vector has to be scaled by -1.0 since the 
                  // load vector for a radiation is assumed +ve for 
                  // flux leaving the surface, which is opposite to the 
                  // FE sign convention.
                  radiation_load_vector.scale(-1.0);
                  
                  DenseVector<double> vector(radiation_load_vector.size());
                  
                  for (unsigned int i=0; i < vector.size(); i++)
                    vector(i) = radiation_load_vector.el(i);
                  
                  vector_it->second->add_vector(vector,
                                                rad_nodal_dofs);
                }
            }
        }
          break;
          
        case FORCE_VECTOR_SENSITIVITY_ENUM_ID:
        {
          if (this->analysis_type_enum_ID == Discipline::NONLINEAR_ANALYSIS::num())
            {
              NumericVector<double>& sol_vec = this->analysis_driver->getCurrentSolution(0);
              FESystemNumerics::PetscSeqVector<double> radiation_load_vector, nodal_temp_vals;
              std::vector<unsigned int> rad_nodal_dofs;
              
              std::vector<RadiationCavityAnalysis*>::iterator it, end;
              it = this->radiation_cavity_analysis_vector.begin();
              end = this->radiation_cavity_analysis_vector.end();
              
              for (; it != end; it++)
                {
                  // get the nodal dofs for the cavity
                  rad_nodal_dofs.clear();
                  (*it)->getFENodalDofIDs(rad_nodal_dofs);
                  
                  // prepare the nodal temperature vector for these dofs
                  nodal_temp_vals.resize(rad_nodal_dofs.size());
                  radiation_load_vector.resize(rad_nodal_dofs.size());
                  
                  for (unsigned int i=0; i < nodal_temp_vals.size(); i++)
                    nodal_temp_vals.set(i, sol_vec(rad_nodal_dofs[i]));
                  
                  (*it)->getFENodalLoadVectorSensitivity
                  (nodal_temp_vals,
                   radiation_load_vector);
                  
                  // the vector has to be scaled by -1.0 since the 
                  // load vector for a radiation is assumed +ve for 
                  // flux leaving the surface, which is opposite to the 
                  // FE sign convention.
                  radiation_load_vector.scale(-1.0);
                  
                  DenseVector<double> vector(radiation_load_vector.size());
                  
                  for (unsigned int i=0; i < vector.size(); i++)
                    vector(i) = radiation_load_vector.el(i);
                  
                  vector_it->second->add_vector(vector,
                                                rad_nodal_dofs);
                }
            }
        }
          break;
          
        default:
        {
          // keep going, nothing to be done here.
        }
      }
    }
}



const NumericVector<double>& 
Discipline::ThermalAnalysis::getSolutionVector(const unsigned int transient_order,
                                               const bool sensitivity)
{
  if (transient_order > 0)
    Assert(this->solution_base->getDisciplineTransientNatureEnumID
           (Discipline::THERMAL_DISCIPLINE::num()) == 
           Discipline::TRANSIENT_SOLUTION::num(), 
           ExcInternalError());
  
  switch (this->analysis_driver->getCurrentAnalysisKind())
  {
    case ANALYSIS_ENUM_ID:
    {
      // asking for sensitivity information does not make sense during analysis 
      Assert(!sensitivity, ExcInternalError());
      
      switch (this->solution_base->getDisciplineTransientNatureEnumID
              (Discipline::THERMAL_DISCIPLINE::num()))
      {
        case STEADY_STATE_SOLUTION_ENUM_ID:
        case TRANSIENT_SOLUTION_ENUM_ID:
        {
          return this->analysis_driver->getCurrentSolution(transient_order);
        }
          break;
          
        default:
          Assert(false, ExcInternalError());
      }
      
    }
      break;
      
    case SENSITIVITY_ANALYSIS_ENUM_ID:
    {
      // asking for sensitivity information does not make sense during sensitivity analysis
      Assert(!sensitivity, ExcInternalError());
      
      // a pointer to the data info object
      static FESystemDatabase::DataInfoBase *data_info_ptr = NULL;
      
      // the solution will have to be loaded for these operations
      switch (this->solution_base->getDisciplineTransientNatureEnumID
              (Discipline::THERMAL_DISCIPLINE::num()))
      {
        case STEADY_STATE_SOLUTION_ENUM_ID:
        {
          static FESystemDatabase::TimeIndependentDataInfo data_info;
          data_info.clear();
          data_info_ptr = &data_info;
        }
          break;
          
        case TRANSIENT_SOLUTION_ENUM_ID:
        {
          static FESystemDatabase::TimeDependentDataInfo data_info;
          data_info.clear();
          Driver::TransientAnalysisDriver* transient_driver = 
          dynamic_cast<Driver::TransientAnalysisDriver*>(this->analysis_driver);
          data_info.setTransientIterationInfo(transient_driver->getCurrentIterationNumber(),
                                              transient_driver->getCurrentAnalysisTime());
          data_info.setOrder(transient_order);
          data_info_ptr = &data_info;
        }
          break;
          
        default:
          Assert(false, ExcInternalError());
      }
      
      data_info_ptr->setDisciplineEnumID(Discipline::THERMAL_DISCIPLINE::num());
      data_info_ptr->setName("Solution");
      data_info_ptr->setLoadCase(this->analysis_driver->getCurrentLoadCase());
      
      return this->loaded_solution.loadSolutionFromDatabase
      (*(this->fesystem_controller.global_data_storage.get()), *data_info_ptr);
    }
      break;
      
    case POST_PROCESS_ENUM_ID:
    default:
      Assert(false, ExcInternalError());
  }
  
}




bool
Discipline::ThermalAnalysis::disciplineHasMatrix(const unsigned int qty_enum_ID) const
{
  bool return_val = false;
  
  switch (qty_enum_ID)
  {
    case TRANSIENT_C1_MATRIX_ENUM_ID:
      return_val = true;
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
Discipline::ThermalAnalysis::getDataInfoForQty(const unsigned int qty_enum_ID)
{
  std::auto_ptr<FESystemDatabase::DataInfoBase> data_info;
  bool set_transient_data = false,
  set_DV_ID = false;
  
  switch (qty_enum_ID)
  {
    case TRANSIENT_C1_MATRIX_ENUM_ID:
    {
      // check the kind of analysis in progress
      switch (this->solution_base->getDisciplineTransientNatureEnumID
              (Discipline::THERMAL_DISCIPLINE::num()))
      {
        case STEADY_STATE_SOLUTION_ENUM_ID:
        {
          // the system matrix is dependent on load case for both linear and 
          // nonlinear analysis, since it uses a load dependent part of convective
          // conductance
          // if the data is for sensitivity, set the
          data_info.reset(new FESystemDatabase::TimeIndependentDataInfo);
          data_info->setName("C_Matrix");
        }
          break;
          
        case TRANSIENT_SOLUTION_ENUM_ID:
        {
          FESystemDatabase::TimeDependentDataInfo *ptr = new FESystemDatabase::TimeDependentDataInfo;
          data_info.reset(ptr);
          data_info->setName("C_Matrix");
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
      
			
    case TRANSIENT_C1_MATRIX_SENSITIVITY_ENUM_ID:
    {
      data_info.reset(this->getDataInfoForQty(Driver::TRANSIENT_C1_MATRIX::num()).release());
      set_DV_ID = true;
    }
      break;
      
    case SYSTEM_MATRIX_ENUM_ID:
    {
      // check the kind of analysis in progress
      switch (this->solution_base->getDisciplineTransientNatureEnumID
              (Discipline::THERMAL_DISCIPLINE::num()))
      {
        case STEADY_STATE_SOLUTION_ENUM_ID:
        {
          // the system matrix is dependent on load case for both linear and 
          // nonlinear analysis, since it uses a load dependent part of convective
          // conductance
          // if the data is for sensitivity, set the
          data_info.reset(new FESystemDatabase::TimeIndependentDataInfo);
          data_info->setName("Kc_plus_Kh_Matrix");
        }
          break;
          
        case TRANSIENT_SOLUTION_ENUM_ID:
        {
          // the system matrix is dependent on load case for both linear and 
          // nonlinear analysis, since it uses a load dependent part of convective
          // conductance
          FESystemDatabase::TimeDependentDataInfo *ptr = new FESystemDatabase::TimeDependentDataInfo;
          data_info.reset(ptr);
          data_info->setName("Kc_plus_Kh_Matrix");
          ptr->setOrder(0);
          set_transient_data = true;
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
      
      
    case FORCE_VECTOR_ENUM_ID:
    {
      // check the kind of analysis in progress
      switch (this->solution_base->getDisciplineTransientNatureEnumID
              (Discipline::THERMAL_DISCIPLINE::num()))
      {
        case STEADY_STATE_SOLUTION_ENUM_ID:
        {
          // the system matrix is dependent on load case for both linear and 
          // nonlinear analysis, since it uses a load dependent part of convective
          // conductance
          data_info.reset(new FESystemDatabase::TimeIndependentDataInfo);
          data_info->setName("F_Vector");
        }
          break;
          
        case TRANSIENT_SOLUTION_ENUM_ID:
        {
          // the system matrix is dependent on load case for both linear and 
          // nonlinear analysis, since it uses a load dependent part of convective
          // conductance
          FESystemDatabase::TimeDependentDataInfo *ptr = new FESystemDatabase::TimeDependentDataInfo;
          data_info.reset(ptr);
          data_info->setName("F_Vector");
          ptr->setOrder(0);
          set_transient_data = true;
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
			
    case JACOBIAN_MATRIX_ENUM_ID:
    {
      // check the kind of analysis in progress
      switch (this->solution_base->getDisciplineTransientNatureEnumID
              (Discipline::THERMAL_DISCIPLINE::num()))
      {
        case STEADY_STATE_SOLUTION_ENUM_ID:
        {
          // the system matrix is dependent on load case for both linear and 
          // nonlinear analysis, since it uses a load dependent part of convective
          // conductance
          data_info.reset(new FESystemDatabase::TimeIndependentDataInfo);
          data_info->setName("Jac_Matrix");
        }
          break;
          
        case TRANSIENT_SOLUTION_ENUM_ID:
        {
          // the system matrix is dependent on load case for both linear and 
          // nonlinear analysis, since it uses a load dependent part of convective
          // conductance
          FESystemDatabase::TimeDependentDataInfo *ptr = new FESystemDatabase::TimeDependentDataInfo;
          data_info.reset(ptr);
          data_info->setName("Jac_Matrix");
          ptr->setOrder(0);
          set_transient_data = true;
        }
          break;
          
        default:
          // at this point of time, no other analysis needs the specification of the elem
          // dof vector
          Assert(false, ExcInternalError());
          
      }
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
      data_info->setDisciplineEnumID(Discipline::THERMAL_DISCIPLINE::num());
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

