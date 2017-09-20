// $Id: SolutionBase.C,v 1.10.6.5 2008-08-21 20:42:16 manav Exp $

// FESystem includes
#include "Solutions/SolutionBase.h"
#include "FESystem/FESystemExceptions.h"
#include "FESystem/FESystemController.h"
#include "FESystem/FESystemNumbers.h"
#include "AnalysisDriver/LinearAnalysisDriver.h"
#include "AnalysisDriver/NonLinearAnalysisDriver.h"
#include "AnalysisDriver/EigenProblemAnalysisDriver.h"
#include "AnalysisDriver/TransientAnalysisDriver.h"
#include "FESystem/AnalysisCase.h"
#include "Properties/PropertyDatabase.h"
#include "DesignData/DesignDatabase.h"  
#include "Discipline/DisciplineInfo.h"
#include "Solvers/LinearSolver.h"
#include "Solvers/NonlinearSolver.h"
#include "Solvers/EigenSolver.h"
#include "Solvers/TransientSolver.h"
#include "Solvers/TransientSolverInfo.h"
#include "Solutions/ThermalSolution.h"
#include "Solutions/TransientThermalSolution.h"
#include "Solutions/LinearStressSolution.h"
#include "Solutions/StructuralVibrationEigenSolution.h"
#include "Solutions/LinearizedBucklingEigenSolution.h"
#include "Solutions/TransientStructuralSolution.h"
#include "Solutions/NonlinearStaticStructuralSolution.h"
#include "Solutions/AeroelasticitySolution.h"
#include "Solvers/LinearSolverInfo.h"
#include "Solvers/NonlinearSolverInfo.h"
#include "Solvers/EigenSolverInfo.h"


Solution::SolutionBase::SolutionBase(const unsigned int sol_enum_ID):
fesystem_controller(NULL),
solution_enum_ID(sol_enum_ID),
ID(FESystemNumbers::InvalidID)
{
  
}


Solution::SolutionBase::~SolutionBase()
{
  
}


unsigned int 
Solution::SolutionBase::getSolutionEnumID() const
{
  return this->solution_enum_ID;
}


std::string 
Solution::SolutionBase::getSolutionEnumName() const
{
  return Solution::SolutionEnum::enumName(this->solution_enum_ID);
}



void 
Solution::SolutionBase::attachFESystemController(FESystem::FESystemController& controller)
{
  this->fesystem_controller = &controller;
}



bool
Solution::SolutionBase::checkParticipatingDiscipline(const unsigned int discipline_enum_id)
{
  // first make sure that the participating disciplines have been initialized
  Assert(this->discipline_enum_IDs.size() > 0,
	 ExcInternalError());

  // now check if the discipline is participating in this solution
  unsigned int n_count = this->discipline_enum_IDs.count(discipline_enum_id);
  
  if (n_count == 0)
    return false;
  else
    return true;

  // should never get here
  Assert(false, ExcInternalError());
  return false;
}



void 
Solution::SolutionBase::addParticipatingDiscipline(const unsigned int discipline_enum_ID)
{
  bool insert_success = 
    this->discipline_enum_IDs.insert(discipline_enum_ID).second;
  Assert(insert_success, ExcInternalError());
}




unsigned int
Solution::SolutionBase::getSolverInfoID(const unsigned int solver_class) const
{
  std::map<unsigned int, unsigned int>::const_iterator it, end;
  it = this->solver_info_ID_map.find(solver_class);
  end = this->solver_info_ID_map.end(); 
  
  Assert(it != end,
         FESystemExceptions::ExcIDDoesNotExist("Solver Class", solver_class));
  
  return it->second;
}



void
Solution::SolutionBase::addSolverID(const unsigned int solver_class,
                                    const unsigned int solver_info_ID)
{
  bool insert_success = 
  this->solver_info_ID_map.insert(std::map<unsigned int, unsigned int>::value_type
                                  (solver_class, solver_info_ID)).second;
  Assert(insert_success,
         FESystemExceptions::ExcDuplicateID("Solver Class", solver_class));
  
}




void
Solution::SolutionBase::linearSolve(const unsigned int discipline_enum)
{

  const Discipline::DisciplineInfo& info = 
  this->fesystem_controller->analysis_case->getAnalysisDisciplineInfo
  (discipline_enum);
  
  // first, make sure that the discipline is for a linear analysis only
  Assert(info.getAnalysisTypeEnumID() == Discipline::LINEAR_ANALYSIS::num(),
         ExcInternalError());
  
  // next, get the discipline and drivers for this analysis
  Discipline::AnalysisDisciplineBase *discipline = 
    this->fesystem_controller->getAnalysisDiscipline(discipline_enum);

  discipline->clear();
  discipline->attachSolution(this);
  
  Driver::AnalysisDriver *driver = 
    this->fesystem_controller->getAnalysisDriver(Driver::LINEAR_ANALYSIS_DRIVER::num());
  
  // attach the discipline and driver to each other
  driver->clear();
  driver->attachAnalysisDiscipline(discipline);
  driver->attachSolution(this);
  
  // get the load cases and design variables to be solved for
  std::auto_ptr<std::vector<DesignData::DesignParameter*> > 
    dv_vector(this->fesystem_controller->design_database->getParameters().release());
  
  // next, peform the analysis
  this->fesystem_controller->initPropertyCardsForGlobalParameters();
  
  // attach a linear solver to the driver and perform the analysis
  unsigned int solver_info_id = this->getSolverInfoID(Solver::LINEAR_SOLVER::num());
  const Solver::LinearSolverInfo& linear_info = dynamic_cast<const Solver::LinearSolverInfo&>
    (this->fesystem_controller->analysis_case->getSolverInfo(solver_info_id));
  std::auto_ptr<Solver::FESystemSolverBase> 
    solver(Solver::createLinearSolver(linear_info).release());

  solver->clear();
  driver->attachSolver(solver.get());
  driver->solveForLoadCases(this->load_case_IDs);
  driver->solveForLoadCaseSensitivity(*(dv_vector.get()), this->load_case_IDs);
  //  driver->postProcess();
  
  // now clear the initializations
  this->fesystem_controller->property_database->clearLocalParameterInitializationForAllCards();
  this->fesystem_controller->property_database->clearGlobalParameterInitializationForAllCards();
  
  // now clear the data structures after use
  solver->clear();
  discipline->clear();
  driver->clear();
}



void
Solution::SolutionBase::nonlinearSolve(const unsigned int discipline_enum)
{  
  const Discipline::DisciplineInfo& info = 
  this->fesystem_controller->analysis_case->getAnalysisDisciplineInfo
  (discipline_enum);

  // first, make sure that the discipline is for a linear analysis only
  Assert(info.getAnalysisTypeEnumID() == Discipline::NONLINEAR_ANALYSIS::num(),
         ExcInternalError());
  
  // next, get the discipline and drivers for this analysis
  Discipline::AnalysisDisciplineBase *discipline = 
    this->fesystem_controller->getAnalysisDiscipline(discipline_enum);

  discipline->clear();
  discipline->attachSolution(this);
  
  Driver::AnalysisDriver *driver = 
    this->fesystem_controller->getAnalysisDriver(Driver::NONLINEAR_ANALYSIS_DRIVER::num());
  
  // attach the discipline and driver to each other
  driver->clear();
  driver->attachAnalysisDiscipline(discipline);
  driver->attachSolution(this);
  
  // get the load cases and design variables to be solved for
  std::auto_ptr<std::vector<DesignData::DesignParameter*> > 
    dv_vector(this->fesystem_controller->design_database->getParameters().release());
  
  // next, peform the analysis
  this->fesystem_controller->initPropertyCardsForGlobalParameters();
  
  // attach a nonlinear solver to the driver and perform the nonlinear analysis
  unsigned int nonlinear_solver_info_id = this->getSolverInfoID(Solver::NONLINEAR_SOLVER::num());
  const Solver::NonlinearSolverInfo& nonlinear_info = 
    dynamic_cast<const Solver::NonlinearSolverInfo&>
    (this->fesystem_controller->analysis_case->getSolverInfo(nonlinear_solver_info_id));

  unsigned int linear_info_ID = nonlinear_info.getLinearSolverInfoID();
  const Solver::LinearSolverInfo& linear_info = 
    dynamic_cast<const Solver::LinearSolverInfo&>
    (this->fesystem_controller->analysis_case->getSolverInfo(linear_info_ID));
  
  std::auto_ptr<Solver::FESystemSolverBase> 
    solver(Solver::createNonlinearSolver(nonlinear_info, linear_info).release());
  
  solver->clear();
  driver->attachSolver(solver.get());
  driver->solveForLoadCases(this->load_case_IDs);
  solver->clear();

  // if there are any DVs, then attach the linear solver and perform 
  // sensitivity analysis
  if (dv_vector->size() > 0)
    {
    linear_info_ID = this->getSolverInfoID(Solver::LINEAR_SOLVER::num());
    const Solver::LinearSolverInfo& sol_linear_info = 
      dynamic_cast<const Solver::LinearSolverInfo&>
      (this->fesystem_controller->analysis_case->getSolverInfo(linear_info_ID));
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
  discipline->clear();
  driver->clear();
}





void
Solution::SolutionBase::eigenSolve(const unsigned int discipline_enum)
{  
  // next, get the discipline and drivers for this analysis
  Discipline::AnalysisDisciplineBase *discipline = 
    this->fesystem_controller->getAnalysisDiscipline(discipline_enum);

  discipline->clear();
  discipline->attachSolution(this);
  
  Driver::AnalysisDriver *driver = 
    this->fesystem_controller->getAnalysisDriver(Driver::EIGENPROBLEM_ANALYSIS_DRIVER::num());
  
  // attach the discipline and driver to each other
  driver->clear();
  driver->attachAnalysisDiscipline(discipline);
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
  discipline->clear();
  driver->clear();
  
}



void
Solution::SolutionBase::linearTransientSolve(const unsigned int discipline_enum)
{  
  const Discipline::DisciplineInfo& info = 
  this->fesystem_controller->analysis_case->getAnalysisDisciplineInfo
  (discipline_enum);
  
  // first, make sure that the discipline is for a linear analysis only
  Assert(info.getAnalysisTypeEnumID() == Discipline::LINEAR_ANALYSIS::num(),
         ExcInternalError());
  
  // next, get the discipline and drivers for this analysis
  Discipline::AnalysisDisciplineBase *discipline = 
    this->fesystem_controller->getAnalysisDiscipline(discipline_enum);
  
  discipline->clear();
  discipline->attachSolution(this);
  
  Driver::AnalysisDriver *driver = this->fesystem_controller->getAnalysisDriver
    (Driver::LINEAR_TRANSIENT_ANALYSIS_DRIVER::num());
  
  // attach the discipline and driver to each other
  driver->clear();
  driver->attachAnalysisDiscipline(discipline);
  driver->attachSolution(this);
  
  // get the load cases and design variables to be solved for
  std::auto_ptr<std::vector<DesignData::DesignParameter*> > 
    dv_vector(this->fesystem_controller->design_database->getParameters().release());
  
  // next, peform the analysis
  this->fesystem_controller->initPropertyCardsForGlobalParameters();
  
  // attach a nonlinear solver to the driver and perform the nonlinear analysis
  unsigned int transient_solver_info_id = this->getSolverInfoID(Solver::LINEAR_TRANSIENT_SOLVER::num());
  const Solver::LinearTransientSolverInfo& transient_info = 
    dynamic_cast<const Solver::LinearTransientSolverInfo&>
    (this->fesystem_controller->analysis_case->getSolverInfo(transient_solver_info_id));
  
  unsigned int linear_info_ID = transient_info.getLinearSolverInfoID();
  const Solver::LinearSolverInfo& linear_info = 
    dynamic_cast<const Solver::LinearSolverInfo&>
    (this->fesystem_controller->analysis_case->getSolverInfo(linear_info_ID));
  
  std::auto_ptr<Solver::FESystemSolverBase> 
    solver(Solver::createLinearTransientSolver(transient_info, linear_info).release());
  
  solver->clear();
  driver->attachSolver(solver.get());
  driver->solveForLoadCases(this->load_case_IDs);

  // if there are any DVs, then attach the linear solver and perform 
  // sensitivity analysis
  if (dv_vector->size() > 0)
    driver->solveForLoadCaseSensitivity(*(dv_vector.get()), this->load_case_IDs);
   
  solver->clear();
  //  driver->postProcess();
  
  // now clear the initializations
  this->fesystem_controller->property_database->clearLocalParameterInitializationForAllCards();
  this->fesystem_controller->property_database->clearGlobalParameterInitializationForAllCards();
  
  // now clear the data structures before returning
  discipline->clear();
  driver->clear();
  
}




void
Solution::SolutionBase::nonlinearTransientSolve(const unsigned int discipline_enum)
{  
  const Discipline::DisciplineInfo& info = 
  this->fesystem_controller->analysis_case->getAnalysisDisciplineInfo
  (discipline_enum);
  
  // first, make sure that the discipline is for a linear analysis only
  Assert(info.getAnalysisTypeEnumID() == Discipline::NONLINEAR_ANALYSIS::num(),
         ExcInternalError());
  
  // next, get the discipline and drivers for this analysis
  Discipline::AnalysisDisciplineBase *discipline = 
    this->fesystem_controller->getAnalysisDiscipline(discipline_enum);
  
  discipline->clear();
  discipline->attachSolution(this);
  
  Driver::AnalysisDriver *driver = this->fesystem_controller->getAnalysisDriver
    (Driver::NONLINEAR_TRANSIENT_ANALYSIS_DRIVER::num());
  
  // attach the discipline and driver to each other
  driver->clear();
  driver->attachAnalysisDiscipline(discipline);
  driver->attachSolution(this);
  
  // get the load cases and design variables to be solved for
  std::auto_ptr<std::vector<DesignData::DesignParameter*> > 
    dv_vector(this->fesystem_controller->design_database->getParameters().release());
  
  // next, peform the analysis
  this->fesystem_controller->initPropertyCardsForGlobalParameters();
  
  // attach a nonlinear solver to the driver and perform the nonlinear analysis
  unsigned int transient_solver_info_id = 
    this->getSolverInfoID(Solver::NONLINEAR_TRANSIENT_SOLVER::num());
  const Solver::NonlinearTransientSolverInfo& transient_info = 
    dynamic_cast<const Solver::NonlinearTransientSolverInfo&>
    (this->fesystem_controller->analysis_case->getSolverInfo(transient_solver_info_id));
  
  unsigned int nonlinear_info_ID = transient_info.getNonlinearSolverInfoID();
  const Solver::NonlinearSolverInfo& nonlinear_info = 
    dynamic_cast<const Solver::NonlinearSolverInfo&>
    (this->fesystem_controller->analysis_case->getSolverInfo(nonlinear_info_ID));

  unsigned int linear_info_ID = nonlinear_info.getLinearSolverInfoID();
  const Solver::LinearSolverInfo& linear_info = 
    dynamic_cast<const Solver::LinearSolverInfo&>
    (this->fesystem_controller->analysis_case->getSolverInfo(linear_info_ID));
  
  std::auto_ptr<Solver::FESystemSolverBase> 
    solver(Solver::createNonlinearTransientSolver(transient_info, nonlinear_info, linear_info).release());
  
  solver->clear();
  driver->attachSolver(solver.get());
  driver->solveForLoadCases(this->load_case_IDs);
  
  // if there are any DVs, then attach the linear solver and perform 
  // sensitivity analysis
  if (dv_vector->size() > 0)
    {
      // get the linear solver info ID
      unsigned int linear_transient_info_ID = 
      this->getSolverInfoID(Solver::LINEAR_TRANSIENT_SOLVER::num());
      
      // get the solver info
      const Solver::LinearTransientSolverInfo& linear_transient_info = 
      dynamic_cast<const Solver::LinearTransientSolverInfo&>
      (this->fesystem_controller->analysis_case->getSolverInfo(linear_transient_info_ID));
      
      // get the linear solver info for this solver
      linear_info_ID = linear_transient_info.getLinearSolverInfoID();
      
      // get the linear solver info 
      const Solver::LinearSolverInfo& linear_info2 = 
      dynamic_cast<const Solver::LinearSolverInfo&>
      (this->fesystem_controller->analysis_case->getSolverInfo(linear_info_ID));
      
      // create the solver object 
      solver.reset(Solver::createLinearTransientSolver(linear_transient_info,
                                                       linear_info2).release());
      solver->clear();

      // get a linear transient analysi driver.
      driver->attachSolver(solver.get());
      driver->solveForLoadCaseSensitivity(*(dv_vector.get()), this->load_case_IDs);
      solver->clear();
    }
  
  solver->clear();
  //  driver->postProcess();
  
  // now clear the initializations
  this->fesystem_controller->property_database->clearLocalParameterInitializationForAllCards();
  this->fesystem_controller->property_database->clearGlobalParameterInitializationForAllCards();
  
  // now clear the data structures before returning
  discipline->clear();
  driver->clear();
  
}







std::auto_ptr<Solution::SolutionBase>
Solution::createSolutionBase(const unsigned int enum_ID)
{
  std::auto_ptr<Solution::SolutionBase> solution(NULL);
  
  switch (enum_ID)
    {
    case THERMAL_SOLUTION_ENUM_ID:
      solution.reset(new Solution::ThermalSolution());
      break;
      
    case TRANSIENT_THERMAL_SOLUTION_ENUM_ID:
      solution.reset(new Solution::TransientThermalSolution());
      break;

    case LINEAR_STRESS_SOLUTION_ENUM_ID:
      solution.reset(new Solution::LinearStressSolution());
      break;

    case STRUCTURAL_VIBRATION_EIGEN_SOLUTION_ENUM_ID:
      solution.reset(new Solution::StructuralVibrationEigenSolution());
      break;

    case LINEARIZED_BUCKLING_EIGEN_SOLUTION_ENUM_ID:
      solution.reset(new Solution::LinearizedBucklingEigenSolution());
      break;

    case TRANSIENT_STRUCTURAL_SOLUTION_ENUM_ID:
      solution.reset(new Solution::TransientStructuralSolution());
      break;

//      case NONLINEAR_STATIC_STRUCTURAL_SOLUTION_ENUM_ID:
//        solution.reset(new Solution::NonlinearStaticStructuralSolution());
//        break;
        
      case AEROELASTICITY_SOLUTION_ENUM_ID:
        solution.reset(new Solution::AeroelasticitySolution());
        break;

      default:
      Assert(false, ExcInternalError());
    }
  
  return solution;
}



