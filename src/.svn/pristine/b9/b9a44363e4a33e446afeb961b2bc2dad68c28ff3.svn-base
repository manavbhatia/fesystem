/*
 *  PistonTheory.cpp
 *  fesystem_xcode
 *
 *  Created by Manav Bhatia on 9/21/08.
 *  Copyright 2008 Manav Bhatia. All rights reserved.
 *
 */

// FESystem include
#include "Discipline/PistonTheory.h"
#include "Discipline/PistonTheoryInfo.h"
#include "AnalysisDriver/TransientAnalysisDriver.h"
#include "AnalysisDriver/AeroelasticityAnalysisDriverBase.h"
#include "Solutions/AeroelasticitySolution.h"
#include "PanelMethods/PistonTheoryTri3.h"
#include "Loads/FlightCondition.h"


Discipline::PistonTheory::PistonTheory(FESystem::FESystemController& controller,
                                       const Discipline::PistonTheoryInfo& info):
Discipline::AerodynamicDisciplineBase(controller, info)
{
  
}


Discipline::PistonTheory::~PistonTheory()
{
  
}


void 
Discipline::PistonTheory::clear()
{
  this->structural_discipline = NULL;
  Discipline::AerodynamicDisciplineBase::clear();
}




bool 
Discipline::PistonTheory::disciplineHasMatrix(const unsigned int qty_enum_ID) const
{
  bool return_val = false;
  
  switch (qty_enum_ID)
  {
    case TRANSIENT_C2_MATRIX_ENUM_ID:
      return_val = false;
      break;
      
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


void
Discipline::PistonTheory::addAnalysisVariables()
{  
  // this discipline does not add any variables of its own, rather 
  // inherits it from the structural analysis. Hence, the dof map from the 
  // structural discipline should be used
  Assert(this->structural_discipline != NULL, ExcInternalError());
  this->dof_map = const_cast<MeshDS::FEDofMap*>(&(this->structural_discipline->getAnalysisDofMap()));
}


void
Discipline::PistonTheory::attachStructuralDiscipline(Discipline::AnalysisDisciplineBase& discipline)
{
  Assert(this->structural_discipline = NULL, ExcInternalError());
  this->structural_discipline = &discipline;
}


unsigned int
Discipline::PistonTheory::getEigenProblemKindEnumID()
{
  Assert(false, ExcInternalError());

  return FESystemNumbers::InvalidID;
}




std::auto_ptr<Loads::DirichletBoundaryConditionDataInfo> 
Discipline::PistonTheory::getBoundaryConditionLoadInfo() const
{  
  // nothign to be done here. 
  Assert(false, ExcInternalError());

  // this discipline has no need for a boundary condition
  std::auto_ptr<Loads::DirichletBoundaryConditionDataInfo> return_ptr;  
  return return_ptr;
}




bool 
Discipline::PistonTheory::calculateTransientMatrix(const unsigned int order,
                                                         bool sensitivity)
{
  // this is to avoid getting warnings about unused parameters
  (void) sensitivity;
  (void) order;

  bool val = false;
  
  switch (order) {
    case 1:
    case 0:
      val = true;
      break;
      
    case 2:
      Assert(false, ExcInternalError());
    default:
      break;
  }
  
  return val;
}




bool 
Discipline::PistonTheory::calculateK(bool sensitivity)
{
  // this is to avoid getting warnings about unused parameters
  (void) sensitivity;
    
  return true;
}





bool 
Discipline::PistonTheory::calculateJac()
{
  // no calculations for a jacobian yet
  Assert(false, ExcInternalError());
  
  // returning a bogus value to avoid warnings of no return
  return false;
}




bool
Discipline::PistonTheory::calculateEigenProblemAMatrix(bool sensitivity) 
{
  // this is to avoid getting warnings about unused parameters
  (void) sensitivity;
  
  // nothing to be done here
  Assert(false, ExcInternalError());
  
  return false;
}




bool
Discipline::PistonTheory::calculateEigenProblemBMatrix(bool sensitivity)
{
  // this is to avoid getting warnings about unused parameters
  (void) sensitivity;

  // nothing to be done here
  Assert(false, ExcInternalError());
  
  return false;
}



bool 
Discipline::PistonTheory::calculateF(bool sensitivity)
{
  // this is to avoid getting warnings about unused parameters
  (void) sensitivity;

  // nothing to be done here
  Assert(false, ExcInternalError());
  
  return false;
}




void 
Discipline::PistonTheory::getElemQty
(FESystemElem::FESystemElemBase* elem,
 std::map<unsigned int, double>& real_elem_data,
 std::map<unsigned int, DenseMatrix<double> >& mat_elem_data,
 std::map<unsigned int, DenseVector<double> >& vec_elem_data,
 const bool sensitivity,
 const unsigned int DV_ID)
{
  // return if the elem type does not belong to this
  if (elem->getFESystemElemTypeEnumID() != PISTON_THEORY_TRI3_ENUM_ID)
    return;
  
  // make sure that valid pointer have been given 
  assert (elem != NULL);
	
  static DenseMatrix<double> matrix;
  static DenseVector<double> vector, elem_dofs, elem_dof_sensitivity;
  static std::vector<DenseVector<double>*> elem_dof_vec(this->getTransientSystemOrder()+1), 
  elem_dof_sens_vec(this->getTransientSystemOrder()+1);
  
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
      surface_load_map[PISTON_THEORY_SURFACE::num()] = SCALAR_SURFACE_LOAD::num();
      created_load_maps = true;
    }
  
  elem_ID = elem->getID();
  load_case = this->analysis_driver->getCurrentLoadCase();
  
  // set the transient nature of this discipline
  static double time_val = 0.0;
  static bool transient_nature = false;
  
  Driver::AeroelasticityAnalysisDriverBase* aeroel_driver = 
  dynamic_cast<Driver::AeroelasticityAnalysisDriverBase*> (this->analysis_driver);
  FESystemElem::PistonTheoryElem* piston_elem = 
  dynamic_cast<FESystemElem::PistonTheoryTri3*>(elem); 
  
  piston_elem->setFluidFlowVector
  (aeroel_driver->getCurrentFlightCondition().getFluidFlowVector());

  
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
  double mach=0.0, sound_speed=0.0;
  double factor1 = 0.0, factor2 = 0.0;
  
  mach = aeroel_driver->getCurrentFlightCondition().getMachNumber(); 
  sound_speed = aeroel_driver->getCurrentFlightCondition().getSpeedOfSound(); 
  
  factor1 = 2.0 / sqrt(mach*mach-1.0) ;
  if (mach > sqrt(2.0))
    factor2 = 2.0 * ((mach*mach-2) / (mach*mach-1)) / sqrt(mach*mach-1.0) / sound_speed;
  else
    factor2 = 2.0 / sqrt(mach*mach-1.0) / sound_speed;
  
  for (; mat_it != mat_end; mat_it++)
    {
      matrix.zero();
      mat_it->second.zero();
      
      switch(mat_it->first)
      {
        case JACOBIAN_MATRIX_ENUM_ID:
        case EIGENPROBLEM_A_MATRIX_ENUM_ID:
        case EIGENPROBLEM_B_MATRIX_ENUM_ID:
        {
          Assert(false, ExcInternalError());
        }
          break;
          
        case SYSTEM_MATRIX_ENUM_ID:
        {
          if (surface_loads.count(PISTON_THEORY_SURFACE::num()) > 0)
            {
              elem->getElementAssembledQty(FESystemElem::PISTON_THEORY_K_MATRIX::num(),
                                           &matrix);
              matrix.scale(factor1);
              mat_it->second += matrix;
            }
        }
          break;
                  
        case SYSTEM_MATRIX_SENSITIVITY_ENUM_ID:
        {
          if (surface_loads.count(PISTON_THEORY_SURFACE::num()) > 0)
            {
              elem->getElementAssembledQtySensitivity(FESystemElem::PISTON_THEORY_K_MATRIX::num(),
                                                      &matrix);
              matrix.scale(factor1);
              mat_it->second += matrix;
            }
        }
          break;                
          
        case TRANSIENT_C1_MATRIX_ENUM_ID:
        {
          if (surface_loads.count(PISTON_THEORY_SURFACE::num()) > 0)
            {
              elem->getElementAssembledQty(FESystemElem::PISTON_THEORY_C1_MATRIX::num(),
                                           &matrix);
              matrix.scale(factor2);
              mat_it->second += matrix;
            }
        }
          break;
          
        case TRANSIENT_C1_MATRIX_SENSITIVITY_ENUM_ID:
        {
          if (surface_loads.count(PISTON_THEORY_SURFACE::num()) > 0)
            {
              elem->getElementAssembledQtySensitivity(FESystemElem::PISTON_THEORY_C1_MATRIX::num(),
                                                      &matrix);
              matrix.scale(factor2);
              mat_it->second += matrix;
            }
        }
          break;                
          

        default:
          abort();
          break;
      }
    }
  
  
  
  // iterate over the vector quantities and get the data
  static std::map<unsigned int, DenseVector<double> >::iterator vec_it, vec_end;
  vec_it = vec_elem_data.begin();
  vec_end = vec_elem_data.end();
  
  for (; vec_it != vec_end; vec_it++)
    {
      vec_it->second.zero();
      
      switch(vec_it->first)
      {
        case FORCE_VECTOR_ENUM_ID:
        case FORCE_VECTOR_SENSITIVITY_ENUM_ID:
        {
          Assert(false, ExcInternalError());
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
        case MODEL_MASS_SENSITIVITY_ENUM_ID:
        default:
          Assert(false, ExcInternalError());
      }
    }      
  
  elem->clearLoadAndDofs();
}




void
Discipline::PistonTheory::performAdditionalOperationsOnGlobalQuantities
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
        case TRANSIENT_C1_MATRIX_ENUM_ID:
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
        case FORCE_VECTOR_SENSITIVITY_ENUM_ID:
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
}



const NumericVector<double>& 
Discipline::PistonTheory::getSolutionVector(const unsigned int transient_order,
                                                  const bool sensitivity)
{
  (void) transient_order;
  (void) sensitivity;
  Assert(false, ExcInternalError());
}




std::auto_ptr<FESystemDatabase::DataInfoBase> 
Discipline::PistonTheory::getDataInfoForQty(const unsigned int qty_enum_ID)
{
  std::auto_ptr<FESystemDatabase::DataInfoBase> data_info;
  bool set_DV_ID = false;
  
  switch (qty_enum_ID)
  {
    case TRANSIENT_C1_MATRIX_ENUM_ID:
    {
      data_info.reset(new FESystemDatabase::TimeIndependentDataInfo);
      data_info->setName("C1_Matrix");
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
      data_info.reset(new FESystemDatabase::TimeIndependentDataInfo);
      data_info->setName("K_Matrix");
    }
      break;
      
      
    case SYSTEM_MATRIX_SENSITIVITY_ENUM_ID:
    {
      data_info.reset(this->getDataInfoForQty(Driver::SYSTEM_MATRIX::num()).release());
      set_DV_ID = true;
    }
      break;
      
      
      
    case TRANSIENT_C2_MATRIX_ENUM_ID:
    case EIGENPROBLEM_A_MATRIX_ENUM_ID:
    case EIGENPROBLEM_B_MATRIX_ENUM_ID:
    case TRANSIENT_C2_MATRIX_SENSITIVITY_ENUM_ID:
    case EIGENPROBLEM_A_MATRIX_SENSITIVITY_ENUM_ID:
    case EIGENPROBLEM_B_MATRIX_SENSITIVITY_ENUM_ID:
    {
      // this discipline does not have an eigensolution defined for it
      Assert(false, ExcInternalError());
    }
      break;
      
      
    case FORCE_VECTOR_ENUM_ID:
    case FORCE_VECTOR_SENSITIVITY_ENUM_ID:
    {
      // this discipline does not have any force defined for it
      Assert(false, ExcInternalError());
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
      data_info->setDisciplineEnumID(Discipline::PISTON_THEORY::num());
      data_info->setLoadCase(this->analysis_driver->getCurrentLoadCase());
    }
  
  
  return data_info;
}




