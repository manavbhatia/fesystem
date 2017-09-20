// $Id: AnalysisDriver.C,v 1.28.6.6 2008/08/21 00:50:52 manav Exp $

// C++ includes



// FESystem inlcudes
#include "AnalysisDriver/AnalysisDriver.h"
#include "FESystem/FESystemController.h"
#include "FESystem/AnalysisCase.h"
#include "Mesh/MeshList.h"
#include "Mesh/FEMeshData.h"
#include "Utilities/Log.h"
#include "Solvers/FESystemSolverBase.h"
#include "DesignData/DesignParameter.h"
#include "Discipline/DisciplineInfo.h"
#include "Loads/LoadDatabase.h"
#include "Loads/LoadCombination.h"
#include "Loads/LoadDataInfo.h"
#include "Utilities/ParallelUtility.h"


// libMesh includes
#include "geom/node.h"
#include "geom/elem.h"
#include "numerics/sparse_matrix.h"
#include "numerics/numeric_vector.h"


Driver::AnalysisDriver::AnalysisDriver(const unsigned int ID,
                                       FESystem::FESystemController& fesys_controller,
                                       const unsigned int driver_type_id):
driver_ID(ID),
fesystem_controller(fesys_controller),
analysis_driver_type_enum_ID(driver_type_id),
solver(NULL),
solution_base(NULL),
current_load_case(0),
current_DV(NULL),
current_analysis_kind_enum_ID(FESystemNumbers::InvalidID)
{

}






Driver::AnalysisDriver::~AnalysisDriver()
{
  this->clear();
}






void
Driver::AnalysisDriver::clear()
{
  // this will need clearing all the data structures:
  this->solver = NULL;

  this->analysis_discipline_map.clear();
  
  this->solution_base = NULL;
  
  this->current_load_case = FESystemNumbers::InvalidID;
  
  this->current_DV = NULL;
  
  this->current_solution.clear();
  
  this->current_analysis_kind_enum_ID = FESystemNumbers::InvalidID;
  
  // -- delete the vectors
  {
    // get the iterators
    std::map<std::string, NumericVector<double>* >::iterator pos = this->vectors.begin();
    std::map<std::string, NumericVector<double>* >::const_iterator end = this->vectors.end();
		
    for (; pos != end; ++pos)
      {
        // if the pointer is not null, clear the data, and delete it
        if (pos->second != NULL)
          {
            pos->second->clear();
            delete pos->second;
            pos->second = NULL;
          }
      }
  }
	
	
  // delete matrices
  {
    // get the iterators
    std::map<std::string, SparseMatrix<double>* >::iterator pos = this->matrices.begin();
    std::map<std::string, SparseMatrix<double>* >::const_iterator end = this->matrices.end();
    
    for (; pos != end; ++pos)
      {
        // if the pointer is not already null, clear the data, and delete it
        if (pos->second != NULL)
          {
            pos->second->clear ();
            delete pos->second;
            pos->second = NULL;
          }
      }
  }
}



void
Driver::AnalysisDriver::attachSolver(Solver::FESystemSolverBase* solver_base)
{
  Assert(solver_base != NULL, ExcInternalError());
  this->solver = solver_base;
  this->solver->attachAnalysisDriver(this);
}




void 
Driver::AnalysisDriver::attachAnalysisDiscipline
(Discipline::AnalysisDisciplineBase* discipline)
{
  Assert(discipline != NULL, ExcEmptyObject());

  // copy the pointer, but make sure that a discipline has not already been attached 
  // for this
  unsigned int enum_ID = discipline->getDisciplineEnumID();

  Assert(this->analysis_discipline_map.count(enum_ID) == 0,
         Driver::AnalysisDriver::ExcInitBeforeClear());
    
  this->analysis_discipline_map.insert(std::map<unsigned int, Discipline::AnalysisDisciplineBase*>::
                                       value_type(enum_ID, discipline));
  
  // attach itself to the discipline
  discipline->attachAnalysisDriver(this);

}



NumericVector<double>& 
Driver::AnalysisDriver::getCurrentSolution(const unsigned int transient_order)
{
  Assert(transient_order < this->current_solution.size(), ExcInternalError());
  Assert(this->current_solution[transient_order] != NULL, ExcEmptyObject());
  return *(this->current_solution[transient_order]);
}



void 
Driver::AnalysisDriver::postProcess(const std::vector<unsigned int>& load_cases, unsigned int disciplin_enum_ID)
{
  // make sure that this discipline exists in the map
  std::map<unsigned int, Discipline::AnalysisDisciplineBase*>::iterator it, end;
  it = this->analysis_discipline_map.find(disciplin_enum_ID);
  end = this->analysis_discipline_map.end();
  
  Assert(it != end, Driver::AnalysisDriver::ExcInitBeforeClear());
  
  // ask the analysis discipline to calculate and store the post
  // process quantity
  it->second->calculatePostProcessQty(load_cases);
}



void
Driver::AnalysisDriver::solveForLoadCases(const std::vector<unsigned int>& load_cases)
{
  // make sure that the solver has been set
  Assert(this->solver != NULL, ExcInternalError());
  
  this->initialize();
  
  // set the analysis type
  this->current_analysis_kind_enum_ID = Driver::ANALYSIS::num();
  
  std::vector<unsigned int>::const_iterator load_case_it = load_cases.begin();
  std::vector<unsigned int>::const_iterator load_case_end = load_cases.end();	
	
  // iterate over each load case, and ask solver to solve for it
  for (; load_case_it != load_case_end; load_case_it++)
    {
      // set current load case
      this->current_load_case = *load_case_it;
      
      // tell solver to solve for the load case
      this->solveCurrentLoadCase();
    }
  
  // now that the solution has been completed, reset the local data
  this->current_load_case = FESystemNumbers::InvalidID;
  this->current_solution.clear();
  this->current_analysis_kind_enum_ID = FESystemNumbers::InvalidID;
}





void 
Driver::AnalysisDriver::solveForLoadCaseSensitivity
(const std::vector<DesignData::DesignParameter*>& design_params,
 const std::vector<unsigned int>& load_cases)
{  
  // if there are no DVs, just return
  if (design_params.size() == 0)
    return;
  
  // make sure that the solver has been set
  Assert(this->solver != NULL, ExcInternalError());
  
  this->initialize();

  // set the analysis kind to sensitivity analysis
  this->current_analysis_kind_enum_ID = Driver::SENSITIVITY_ANALYSIS::num();	
  
  std::vector<unsigned int>::const_iterator load_case_it = load_cases.begin();
  std::vector<unsigned int>::const_iterator load_case_end = load_cases.end();	
  
  std::vector<DesignData::DesignParameter*>::const_iterator dv_it, dv_end;
  dv_end = design_params.end();
	
  
  // iterate over each load case, and ask solver to solve for it
  for (; load_case_it != load_case_end; load_case_it++)
    {
      dv_it = design_params.begin();
      
      for (; dv_it != dv_end; dv_it++)
        {
          // set current load case
          this->current_load_case = *load_case_it;
          this->current_DV = *dv_it;
          
          // solve for this sensitivity analysis
          this->solveCurrentLoadCase();
        }
    }
  
  // now that the solution has been completed, reset the local data
  this->current_load_case = FESystemNumbers::InvalidID;
  this->current_solution.clear();
  this->current_analysis_kind_enum_ID = FESystemNumbers::InvalidID;
  
  this->current_DV = NULL;
}






NumericVector<double> & 
Driver::AnalysisDriver::addVector(const std::string& vec_name)
{
  // Return the vector if it is already there.
  if (this->haveVector(vec_name))
    {
      return *(this->vectors[vec_name]);
    }
	
  // Otherwise build the vector and return it.
  NumericVector<double>* buf = NumericVector<double>::build().release();
  this->vectors.insert (std::make_pair (vec_name, buf));
	
  return *buf;
}






SparseMatrix<double> & 
Driver::AnalysisDriver::addMatrix(const std::string& mat_name)
{
  // Return the matrix if it is already there.
  if (this->haveMatrix(mat_name))
    return *(this->matrices[mat_name]);
	
  // Otherwise build the matrix and return it.
  SparseMatrix<double>* buf = SparseMatrix<double>::build().release();
  this->matrices.insert (std::make_pair (mat_name, buf));
	
  return *buf;
}








void
Driver::AnalysisDriver::writeSolutionVectorToLog
( FESystemDatabase::GenericDataInfoBase& data_info,
 const NumericVector<double>& sol, unsigned int discipline_enum_ID)
{
  
  // make sure that this discipline exists in the map
  std::map<unsigned int, Discipline::AnalysisDisciplineBase*>::iterator discipline_it, discipline_end;
  discipline_it = this->analysis_discipline_map.find(discipline_enum_ID);
  discipline_end = this->analysis_discipline_map.end();
  
  Assert(discipline_it != discipline_end, Driver::AnalysisDriver::ExcInitBeforeClear());
  
  const MeshDS::FEMesh& mesh = discipline_it->second->getAnalysisMesh();
  const MeshDS::FEMeshData& mesh_data = discipline_it->second->getAnalysisMeshData();
    
  if (FESystem::local_processor == 0)
    this->fesystem_controller.log_file->write2(data_info);

  std::vector<double> localized_vector;
  sol.localize_to_one(localized_vector, 0);
    
  if (FESystem::local_processor == 0)
    {
      // now write the dofs to the file
      MeshBase::const_node_iterator           node  = mesh.nodes_begin();
      const MeshBase::const_node_iterator end_node  = mesh.nodes_end();
      
      unsigned int node_id = 0, dof_number=0;
      unsigned int n_vars = discipline_it->second->getNVars();
      std::ostringstream values;
      
      for ( ; node != end_node ; ++node)
        {
          node_id = mesh_data.getInternalIDFromNode(*node);

          values << node_id << "  " ;
          for (unsigned int i=0;i<3; i++)
            values << std::setprecision(16) << std::showpoint << (**node)(i) << "  "; 
          
          for (unsigned int var_it =0; var_it < n_vars; var_it++)
            {
              dof_number = (*node)->dof_number(0,var_it,0);
              values << std::showpoint << localized_vector[dof_number] << "  ";
            }
          values << std::endl;
        }

      const std::string& val_str = values.str();
      this->fesystem_controller.log_file->write(val_str);
    }

}



void
Driver::AnalysisDriver::writeComplexSolutionVectorToLog
( FESystemDatabase::GenericDataInfoBase& data_info, const NumericVector<double>& real_sol,
 const NumericVector<double>& img_sol, unsigned int discipline_enum_ID)
{
  
  // make sure that this discipline exists in the map
  std::map<unsigned int, Discipline::AnalysisDisciplineBase*>::iterator discipline_it, discipline_end;
  discipline_it = this->analysis_discipline_map.find(discipline_enum_ID);
  discipline_end = this->analysis_discipline_map.end();
  
  Assert(discipline_it != discipline_end, Driver::AnalysisDriver::ExcInitBeforeClear());
  
  const MeshDS::FEMesh& mesh = discipline_it->second->getAnalysisMesh();
  const MeshDS::FEMeshData& mesh_data = discipline_it->second->getAnalysisMeshData();
  
  if (FESystem::local_processor == 0)
    this->fesystem_controller.log_file->write2(data_info);
  
  std::vector<double> localized_vector_real;
  real_sol.localize_to_one(localized_vector_real, 0);

  std::vector<double> localized_vector_img;
  img_sol.localize_to_one(localized_vector_img, 0);

  if (FESystem::local_processor == 0)
    {
      // now write the dofs to the file
      MeshBase::const_node_iterator           node  = mesh.nodes_begin();
      const MeshBase::const_node_iterator end_node  = mesh.nodes_end();
      
      unsigned int node_id = 0, dof_number=0;
      unsigned int n_vars = discipline_it->second->getNVars();
      std::ostringstream values;
      
      for ( ; node != end_node ; ++node)
        {
          node_id = mesh_data.getInternalIDFromNode(*node);
          
          values << node_id << "  " ;
          for (unsigned int i=0;i<3; i++)
            values << std::setprecision(16) << std::showpoint << (**node)(i) << "  "; 
          
          for (unsigned int var_it =0; var_it < n_vars; var_it++)
            {
              dof_number = (*node)->dof_number(0,var_it,0);
              values << std::showpoint 
              << localized_vector_real[dof_number] << "  + i "
              << localized_vector_img[dof_number] ;
            }
          values << std::endl;
        }
      
      const std::string& val_str = values.str();
      this->fesystem_controller.log_file->write(val_str);
    }
  
}




void 
Driver::AnalysisDriver::applyBoundaryConditionsToVectorForCurrentAnalysis
(NumericVector<double>& vec, 
 const bool if_apply_bc_vals, unsigned int discipline_enum_ID)
{
  // make sure that this discipline exists in the map
  std::map<unsigned int, Discipline::AnalysisDisciplineBase*>::iterator it, end;
  it = this->analysis_discipline_map.find(discipline_enum_ID);
  end = this->analysis_discipline_map.end();
  
  Assert(it != end, Driver::AnalysisDriver::ExcInitBeforeClear());

  LoadDatabase& load_database = *(this->fesystem_controller.load_database.get());
  const MeshDS::FEMeshData& mesh_data = it->second->getAnalysisMeshData();
  
  // get the boundary conditions
  std::auto_ptr<Loads::DirichletBoundaryConditionDataInfo> 
  load_info(it->second->getBoundaryConditionLoadInfo().release());
  
  FESystemUtility::AutoPtrVector<Loads::DirichletBoundaryConditionCombination> 
  bcond(load_database.getAllLoadCombinations
        <Loads::DirichletBoundaryConditionDataInfo, Loads::DirichletBoundaryConditionCombination>
        (*load_info).release());
  
  // to apply the boundary condition, get the global dof_number for each
  // constrained dof and use the penalty method to apply the boundary condition
  std::vector<Loads::DirichletBoundaryConditionCombination*>::const_iterator bc_it, bc_end;
  bc_it = bcond.get()->begin();
  bc_end = bcond.get()->end();
  
  unsigned int constrained_node, constrained_dof, dof;
  double value = 0.0;
  
  for (; bc_it != bc_end ; bc_it++)
    {
      constrained_node =  (**bc_it).getNodeID();
      constrained_dof = (**bc_it).getDofNumber();
      
      dof = 
      (mesh_data.getNodeFromForeignID(constrained_node))->dof_number(0,constrained_dof-1,0);
      
      if (if_apply_bc_vals)
        {
          value =  (**bc_it).getValue();
          vec.set(dof, value);	    
        }
      else
        vec.set(dof, 0.0);
    } 
}




void 
Driver::AnalysisDriver::applyBoundaryConditionsToMatrixForCurrentAnalysis
(SparseMatrix<double>& mat,
 double value,
 unsigned int discipline_enum_ID)
{
  // make sure that this discipline exists in the map
  std::map<unsigned int, Discipline::AnalysisDisciplineBase*>::iterator it, end;
  it = this->analysis_discipline_map.find(discipline_enum_ID);
  end = this->analysis_discipline_map.end();
  
  Assert(it != end, Driver::AnalysisDriver::ExcInitBeforeClear());

  LoadDatabase& load_database = *(this->fesystem_controller.load_database.get());
  const MeshDS::FEMeshData& mesh_data = it->second->getAnalysisMeshData();
  
  // get the boundary conditions
  std::auto_ptr<Loads::DirichletBoundaryConditionDataInfo> 
  load_info(it->second->getBoundaryConditionLoadInfo().release());
  
  FESystemUtility::AutoPtrVector<Loads::DirichletBoundaryConditionCombination> 
  bcond(load_database.getAllLoadCombinations
        <Loads::DirichletBoundaryConditionDataInfo, Loads::DirichletBoundaryConditionCombination>
        (*load_info).release());
  
  // to apply the boundary condition, get the global dof_number for each
  // constrained dof and use the penalty method to apply the boundary condition
  std::vector<Loads::DirichletBoundaryConditionCombination*>::const_iterator bc_it, bc_end;
  bc_it = bcond.get()->begin();
  bc_end = bcond.get()->end();
  
  std::vector<unsigned int> bc_rows;
  unsigned int constrained_node, constrained_dof, dof;
  
  for (; bc_it != bc_end ; bc_it++)
    {
      constrained_node =  (**bc_it).getNodeID();
      constrained_dof = (**bc_it).getDofNumber();
      dof = 
      (mesh_data.getNodeFromForeignID(constrained_node))->dof_number(0,(constrained_dof-1),0);
	    
      //matrix.add(dof,dof,penalty);
      bc_rows.push_back(dof);
    } 
  
  mat.zero_rows(bc_rows, value);	
}




void
Driver::AnalysisDriver::localizeSolution(const DofMap& dof_map,
                                         NumericVector<double>& source_vec, 
                                         NumericVector<double>& localized_vec)
{
  // make sure that the source vec and the localized vec have the correct dimensions
  Assert(source_vec.size() == dof_map.n_dofs(), ExcInternalError());
  Assert(source_vec.size() == localized_vec.size(), ExcInternalError());
  
  // get the send list from the dofmap
  const std::vector<unsigned int>& send_list = dof_map.get_send_list();
  
  source_vec.close();
  localized_vec.zero();
  localized_vec.close();
  
  source_vec.localize(localized_vec, send_list);
}



