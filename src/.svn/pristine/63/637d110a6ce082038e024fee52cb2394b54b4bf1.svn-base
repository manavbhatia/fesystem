/*
 *  AeroelasticitySolution.C
 *  fesystem_xcode
 *
 *  Created by Manav Bhatia on 9/21/08.
 *  Copyright 2008 Manav Bhatia. All rights reserved.
 *
 */

#include "AeroelasticitySolution.h"


// $Id: AeroelasticitySolution.C,v 1.5.6.4 2008/08/25 04:39:11 manav Exp $

// FESystem includes
#include "Solutions/AeroelasticitySolution.h"
#include "Utilities/InputOutputUtility.h"
#include "Discipline/StructuralAnalysis.h"
#include "Discipline/ThermalAnalysis.h"
#include "AnalysisDriver/AnalysisDriver.h"
#include "DesignData/DesignDatabase.h"
#include "Solvers/EigenSolver.h"
#include "Solvers/EigenSolverInfo.h"
#include "Solvers/LinearSolver.h"
#include "Solvers/LinearSolverInfo.h"
#include "Database/DataInfo.h"
#include "FESystem/AnalysisCase.h"
#include "Properties/PropertyDatabase.h"
#include "FESystem/FESystemExceptions.h"
#include "Loads/FlightCondition.h"


Solution::AeroelasticitySolution::AeroelasticitySolution():
Solution::SolutionBase(Solution::AEROELASTICITY_SOLUTION::num()),
use_modal_order_reduction(false),
coupled_thermoelastic(false),
include_geometric_effects(false),
completed_linear_solve(false),
interpolation_case_ID(0),
n_lag_terms(0),
aerodynamic_discipline_enum_ID(0),
aeroelastic_driver_enum_ID(0)
{
  // add the participating discipline to this solution
  this->addParticipatingDiscipline(Discipline::STRUCTURAL_DISCIPLINE::num());
}



Solution::AeroelasticitySolution::~AeroelasticitySolution()
{
  
}




unsigned int 
Solution::AeroelasticitySolution::getDisciplineTransientNatureEnumID(const unsigned int discipline) const 
{
  unsigned int enum_ID = 0;
  
  switch (discipline)
  {
    case STRUCTURAL_DISCIPLINE_ENUM_ID:
      enum_ID = Discipline::STEADY_STATE_SOLUTION::num();
      break;
      
    default:
      Assert(false, ExcInternalError());
  }
  
  return enum_ID;
}



unsigned int
Solution::AeroelasticitySolution::numberOfAerodynamicLagTerms() const
{
  return this->n_lag_terms;
}



unsigned int
Solution::AeroelasticitySolution::getAeroelasticDriverEnumID() const
{
  return this->aeroelastic_driver_enum_ID;
}




unsigned int
Solution::AeroelasticitySolution::getAerodynamicDisciplineEnumID() const
{
  return this->aerodynamic_discipline_enum_ID;
}



const std::vector<Loads::FlightCondition>&
Solution::AeroelasticitySolution::getFlightConditionVector() const
{
  return this->flight_condition_vector;
}


bool 
Solution::AeroelasticitySolution::ifIncludeGeometricEffects() const
{
  bool val = false;
  
  if (this->include_geometric_effects && this->completed_linear_solve)
    val = true;

  return val;
}


void
Solution::AeroelasticitySolution::solve()
{
  // solve the linear system if the geometric effects are needed
  if (this->include_geometric_effects)
    this->linearSolve(Discipline::STRUCTURAL_DISCIPLINE::num());
  
  this->completed_linear_solve = true;
  
  // solve the eigen problem
  this->aeroelasticSolve();
}


void
Solution::AeroelasticitySolution::aeroelasticSolve()
{  
  // next, get the discipline and drivers for this analysis
  Discipline::AnalysisDisciplineBase * structural_discipline = 
  this->fesystem_controller->getAnalysisDiscipline(Discipline::STRUCTURAL_DISCIPLINE::num());
  structural_discipline->clear();
  structural_discipline->attachSolution(this);  

  Discipline::AnalysisDisciplineBase * aero_discipline = 
  this->fesystem_controller->getAnalysisDiscipline(this->getAerodynamicDisciplineEnumID());
  aero_discipline->clear();
  aero_discipline->attachSolution(this);  
  
  // get the driver and initialize it
  unsigned int driver_enum_ID = this->getAeroelasticDriverEnumID();
  
  Driver::AnalysisDriver *driver = 
  this->fesystem_controller->getAnalysisDriver(driver_enum_ID);
  
  // attach the discipline and driver to each other
  driver->clear();
  driver->attachAnalysisDiscipline(structural_discipline);
  driver->attachAnalysisDiscipline(aero_discipline);
  driver->attachSolution(this);
  
  // get the load cases and design variables to be solved for
  std::auto_ptr<std::vector<DesignData::DesignParameter*> > 
  dv_vector(this->fesystem_controller->design_database->getParameters().release());
  
  // next, peform the analysis
  this->fesystem_controller->initPropertyCardsForGlobalParameters();
  
  // attach a nonlinear solver to the driver and perform the nonlinear analysis
  unsigned int eigen_solver_info_id = this->getSolverInfoID(Solver::EIGEN_SOLVER::num());
  const Solver::EigenSolverInfo& eigen_info = 
  dynamic_cast<const Solver::EigenSolverInfo&>
  (this->fesystem_controller->analysis_case->getSolverInfo(eigen_solver_info_id));
  
  unsigned int linear_info_ID = eigen_info.getLinearSolverInfoID();
  const Solver::LinearSolverInfo& linear_info = 
  dynamic_cast<const Solver::LinearSolverInfo&>
  (this->fesystem_controller->analysis_case->getSolverInfo(linear_info_ID));
  
  std::auto_ptr<Solver::FESystemSolverBase>
  solver(Solver::createEigenSolver(eigen_info, linear_info).release());
  
  solver->clear();
  
  driver->attachSolver(solver.get());
  driver->solveForLoadCases(this->load_case_IDs);
  solver->clear();
  
  // if there are any DVs, then attach the linear solver and perform 
  // sensitivity analysis
  if (dv_vector->size() > 0)
    {
      unsigned int solver_info_id = this->getSolverInfoID(Solver::LINEAR_SOLVER::num());
      const Solver::LinearSolverInfo& sol_linear_info =
      dynamic_cast<const Solver::LinearSolverInfo&>
      (this->fesystem_controller->analysis_case->getSolverInfo(solver_info_id));
      solver.reset(Solver::createLinearSolver(sol_linear_info).release());
      
      solver->clear();
      driver->attachSolver(solver.get());
      driver->solveForLoadCaseSensitivity(*(dv_vector.get()), this->load_case_IDs);
      solver->clear();
    }
  //  driver->postProcess();
  
  // now clear the initializations
  this->fesystem_controller->property_database->clearLocalParameterInitializationForAllCards();
  this->fesystem_controller->property_database->clearGlobalParameterInitializationForAllCards();
  
  // now clear the data structures before returning
  structural_discipline->clear();
  aero_discipline->clear();
  driver->clear();
  
}



std::istream&
Solution::AeroelasticitySolution::readFromInputStream(std::istream& input)
{
  std::string tag;
  unsigned int num = 0;
  
  FESystemIOUtility::readFromInput(input, Solution::AEROELASTICITY_SOLUTION::name());
  FESystemIOUtility::readFromInput(input, "BEGIN");
  
  FESystemIOUtility::readFromInput(input, "ID", this->ID);

  // read the aerodynamic discipline enum 
  tag.clear();
  FESystemIOUtility::readFromInput(input, "AERODYNAMIC_DISCIPLINE", tag);
  this->aerodynamic_discipline_enum_ID = Discipline::AnalysisDisciplineEnum::enumID(tag);

  // read the aeroelasticity solution kind enum 
  tag.clear();
  FESystemIOUtility::readFromInput(input, "SOLUTION_METHOD", tag);
  this->aeroelastic_driver_enum_ID = Driver::AnalysisDriverTypeEnum::enumID(tag);
  
  
  FESystemIOUtility::readFromInput(input, "COUPLED_THERMAL_STRUCTURAL", this->coupled_thermoelastic);

  FESystemIOUtility::readFromInput(input, "GEOMETRIC_EFFECTS", this->include_geometric_effects);
  
  if (this->coupled_thermoelastic)
    {
      this->addParticipatingDiscipline(Discipline::THERMAL_DISCIPLINE::num());
      FESystemIOUtility::readFromInput(input, "INTERPOLATION_CASE_ID", 
                                       this->interpolation_case_ID);
    }
  
  // Read in the list of load cases
  unsigned int case_ID =0;
  FESystemIOUtility::readFromInput(input, "N_ANALYSIS_LOAD_CASES", num);
  for (unsigned int case_it =0; case_it < num; case_it++)
    {
      case_ID = 0;
      input >> case_ID;
      
      this->load_case_IDs.push_back(case_ID);
    }
  
  // read the dynamic pressure and mach numbers
  double val = 0.0;
  Loads::FlightCondition flt_cond;
  Point vec;
  
  FESystemIOUtility::readFromInput(input, "N_FLIGHT_CONDITIONS", num);
  for (unsigned int i=0; i < num; i++)
    {
      FESystemIOUtility::readFromInput(input, "MACH_NUMBER", val);
      flt_cond.setMachNumber(val);
      FESystemIOUtility::readFromInput(input, "DENSITY", val);
      flt_cond.setDensity(val);
      FESystemIOUtility::readFromInput(input, "DYNAMIC_PRESSURE", val);
      flt_cond.setDynamicPressure(val);
      FESystemIOUtility::readFromInput(input, "FLUID_FLOW_VECTOR");
      for (unsigned int i=0; i < 3; i++)
        {
          input >> val;
          vec(i) = val;
        }
      flt_cond.setFluidFlowVector(vec);
      this->flight_condition_vector.push_back(flt_cond);
    }
  
  
  
  FESystemIOUtility::readFromInput(input, "N_SOLVER_INFO_ID", num);
  unsigned int solver_class = 0, info_id = 0;
  for (unsigned int i=0; i < num; i++)
    {
      tag.clear();
      input >> tag;
      solver_class = Solver::SolverClassEnum::enumID(tag);
      input >> info_id;
      
      this->addSolverID(solver_class, info_id);
    }
  
  FESystemIOUtility::readFromInput(input, Solution::AEROELASTICITY_SOLUTION::name());
  FESystemIOUtility::readFromInput(input, "END");
  
  return input;
}



FESystemUtility::AutoPtrVector<FESystemDatabase::DataInfoBase>
Solution::AeroelasticitySolution::getTimeIndependentNodalSolutionDataInfo
(const unsigned int discipline_enum_ID)
{
  FESystemUtility::AutoPtrVector<FESystemDatabase::DataInfoBase> data_info_vec;
  
  switch (discipline_enum_ID)
  {
    case STRUCTURAL_DISCIPLINE_ENUM_ID:
    {
      std::auto_ptr<std::vector<DesignData::DesignParameter*> > 
      dv_vector(this->fesystem_controller->design_database->getParameters().release());
      
      FESystemDatabase::TimeIndependentDataInfo eig_val_data_info;
      
      FESystemDatabase::DataInfoBase *data_info  = NULL;
      
      DenseVector<double> eig_value_real;

      const std::vector<Loads::FlightCondition>& flt_conds = this->getFlightConditionVector();

      for (unsigned int ii=0; ii<flt_conds.size(); ii++)
        {
          std::ostringstream oss_eigval_name_real;
          oss_eigval_name_real << "RogerEigenValue_Real_FltCond_" << ii;
          
          // first get the number of converged eigen values for this load case
          eig_val_data_info.clear();
          eig_val_data_info.setDisciplineEnumID(Discipline::STRUCTURAL_DISCIPLINE::num());
          eig_val_data_info.setName(oss_eigval_name_real.str());
          eig_val_data_info.setLoadCase(1);
          this->fesystem_controller->global_data_storage->fillVector
          (eig_val_data_info, eig_value_real);
          
          for (unsigned int j=0; j < eig_value_real.size(); j++)
            {
              std::ostringstream oss_str_eigvec_name_real;
              oss_str_eigvec_name_real << "EigenVector_Real_FltCond_" << ii;
              
              data_info = new FESystemDatabase::TimeIndependentEigenDataInfo;
              
              data_info->setDisciplineEnumID(Discipline::STRUCTURAL_DISCIPLINE::num());
              data_info->setName(oss_str_eigvec_name_real.str());
              data_info->setLoadCase(1);
              dynamic_cast<FESystemDatabase::TimeIndependentEigenDataInfo*>(data_info)->setModeInfo(j+1);
              
              data_info_vec.push_back(data_info);
              data_info = NULL;
            }
          
//          for (unsigned int i=0; i < this->load_case_IDs.size(); i++)
//            for (unsigned int j=0; j < dv_vector->size(); j++)
//              {
//                // first get the number of converged eigen values for this load case
//                eig_val_data_info.clear();
//                eig_val_data_info.setDisciplineEnumID(Discipline::STRUCTURAL_DISCIPLINE::num());
//                eig_val_data_info.setName("EigenValuesReal");
//                eig_val_data_info.setLoadCase(this->load_case_IDs[i]);
//                eig_val_data_info.setDVID((*dv_vector)[j]->getID());
//                this->fesystem_controller->global_data_storage->fillVector(eig_val_data_info,
//                                                                           eig_value_real);
//                
//                for (unsigned int k=0; k < eig_value_real.size(); k++)
//                  {
//                    data_info = new FESystemDatabase::TimeIndependentEigenDataInfo;
//                    
//                    data_info->setDisciplineEnumID(Discipline::STRUCTURAL_DISCIPLINE::num());
//                    data_info->setName("EigenVector_Real");
//                    data_info->setLoadCase(this->load_case_IDs[i]);
//                    data_info->setDVID((*dv_vector)[j]->getID());
//                    dynamic_cast<FESystemDatabase::TimeIndependentEigenDataInfo*>(data_info)->setModeInfo(k+1);
//                    
//                    data_info_vec.push_back(data_info);
//                    data_info = NULL;
//                  }
//              }
        }
    }
      break;
      
    default:
      Assert(false, ExcInternalError());
  }
  
  return data_info_vec;
}


FESystemUtility::AutoPtrVector<FESystemDatabase::TransientDataInfoSet>
Solution::AeroelasticitySolution::getTimeDependentNodalSolutionDataInfo
(const unsigned int discipline_enum_ID)
{
  // unused parameter
  (void) discipline_enum_ID;
  
  // this solution does not have any transients. Hence, return an empty vector
  FESystemUtility::AutoPtrVector<FESystemDatabase::TransientDataInfoSet> vec;
  return vec;
}





