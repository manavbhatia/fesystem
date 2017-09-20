// $Id: DirectInterpolation.C,v 1.1.2.1 2007-06-13 14:57:20 manav Exp $

// C++ includes


// FESystem includes
#include "Interpolation/DirectInterpolation.h"
#include "Interpolation/InterpolationCase.h"
#include "FESystem/FESystemController.h"
#include "FESystem/AnalysisCase.h"
#include "FESystem/FESystemExceptions.h"
#include "Database/GlobalDataStorage.h"
#include "DesignData/DesignParameter.h"
#include "DesignData/DesignDatabase.h"
#include "Solutions/SolutionBase.h"
#include "Discipline/AnalysisDisciplineBase.h"
#include "Discipline/ThermalAnalysis.h"
#include "Discipline/StructuralAnalysis.h"
#include "Discipline/DisciplineInfo.h"
#include "Numerics/PetscSeqVector.h"
#include "Loads/LoadDatabase.h"
#include "Loads/LoadDataInfo.h"
#include "Loads/LoadSet.h"
#include "Loads/LoadCase.h"
#include "Utilities/InputOutputUtility.h"

// libMesh includes


// Forward declerations



DirectInterpolation::DirectInterpolation(FESystem::FESystemController& controller,
                                         const InterpolationCase& interp_case):
InterpolationBase(controller, interp_case)
{
  
}



DirectInterpolation::~DirectInterpolation()
{
  
}




void
DirectInterpolation::interpolate()
{
  this->createLoads();
}




void
DirectInterpolation::createLoads()
{
  unsigned int load_ID = 0, node_ID = 0, dof_number = 0;
  
  
  // reference to load database and global database
  LoadDatabase& load_database = 
    *(this->fesystem_controller.load_database.get());
  MeshDS::MeshList& mesh_list = 
    *(this->fesystem_controller.mesh_list.get());
  
  std::auto_ptr<NumericVector<double> > 
    solution(NumericVector<double>::build().release());

  
  // get the load cases and the design variables for the analysis case. 
  std::auto_ptr<std::vector<DesignData::DesignParameter*> > 
    dv_vector(this->fesystem_controller.design_database->getParameters().release());
  
  std::vector<Solution::SolutionBase*> solutions = 
    this->fesystem_controller.analysis_case->getSolutions();
  
  // at this point of time, only a single solution is being supported
  Assert(solutions.size() == 1, ExcInternalError());
  
  // now get the load case
  const std::vector<unsigned int>& load_cases = solutions[0]->getLoadCaseIDs();
  
  // get the discipline info for the involved disciplines in this interpolation
  // for now, it is thermal and structural
  // make sure that the two disciplines are a part of the solution
  unsigned int from_discipline_enum_ID = this->interpolation_case.fromDiscipline(),
    to_discipline_enum_ID = this->interpolation_case.toDiscipline();
    
  Assert(from_discipline_enum_ID == Discipline::THERMAL_DISCIPLINE::num(),
         ExcInternalError());
  Assert(to_discipline_enum_ID == Discipline::STRUCTURAL_DISCIPLINE::num(),
         ExcInternalError());
  Assert(solutions[0]->checkParticipatingDiscipline(from_discipline_enum_ID),
         ExcInternalError());
  Assert(solutions[0]->checkParticipatingDiscipline(to_discipline_enum_ID),
         ExcInternalError());
  
  // get the discipline infos
  const Discipline::DisciplineInfo
    &from_discipline_info = 
    this->fesystem_controller.analysis_case->getAnalysisDisciplineInfo(from_discipline_enum_ID),
    &to_discipline_info = 
    this->fesystem_controller.analysis_case->getAnalysisDisciplineInfo(to_discipline_enum_ID);
  
  // get the mesh for the from and to discipline
  const MeshDS::FEMesh
    *from_mesh = mesh_list.getMeshFromID(from_discipline_info.getMeshID()),
    *to_mesh = mesh_list.getMeshFromID(to_discipline_info.getMeshID());
  const MeshDS::FEDofMap
    *from_dof_map = mesh_list.getDofMapFromID(from_discipline_info.getMeshID());
  const MeshDS::FEMeshData* 
    from_mesh_data = mesh_list.getMeshDataFromID(from_discipline_info.getMeshID());
  
  // check the analysis type for the "from" discipline. If it is a static analysis, then
  // only one solution needs to be imported. For transient analysis, all stored time values will 
  // have to be imported. 
  unsigned int from_discipline_transient_nature = 
    solutions[0]->getDisciplineTransientNatureEnumID(from_discipline_enum_ID),
    to_discipline_transient_nature = 
    solutions[0]->getDisciplineTransientNatureEnumID(to_discipline_enum_ID);
  

  // create a data info here to get the solution 
  if (from_discipline_transient_nature == Discipline::TRANSIENT_SOLUTION::num())
    {
    std::auto_ptr<FESystemDatabase::TimeDependentDataInfo> data_info
    (new FESystemDatabase::TimeDependentDataInfo);
    
    for (unsigned int i=0; i < load_cases.size(); i++) // iterate over the load cases
      {
      for (unsigned int k=0; k < (dv_vector->size() + 1); k++) // iterate over the load cases
        {
        FESystemNumerics::PetscSeqVector<double> time_vec;
        FESystemDatabase::TimeIndependentDataInfo time_vals_data_info;
        time_vals_data_info.setDisciplineEnumID(from_discipline_enum_ID);
        time_vals_data_info.setName("TimeValues");
        time_vals_data_info.setLoadCase(load_cases[i]);
        this->fesystem_controller.global_data_storage->fillVector(time_vals_data_info, time_vec);
        
        // iterate over the time instants
        for (unsigned int j=0; j < time_vec.size(); j++)
          {
          // clear the data info and set the data for the solution vector 
          // to be read in
          data_info->clear();
          data_info->setDisciplineEnumID(from_discipline_enum_ID);
          data_info->setName("Solution");
          data_info->setOrder(0); // only the zeroth order will be written as load.
          if (k>0) // write the dv value
            data_info->setDVID((*dv_vector)[k-1]->getID());
          data_info->setLoadCase(load_cases[i]);
          data_info->setTransientIterationInfo(j, time_vec.el(j));
          
          this->fesystem_controller.global_data_storage->fillVector(*data_info, *solution);
          
          // make sure that the number of dofs in the mesh are the same as the number of dofs 
          // in this solution vector
          Assert(solution->size() == from_dof_map->n_dofs(), ExcInternalError());


          // iterate over all the dofs, and write the load values to the output stream
          // create a load data info for the loads
          const TimeDependentLoadCase& load_case = dynamic_cast<const TimeDependentLoadCase&>
            (load_database.getLoadCaseFromID(load_cases[i]));

          std::auto_ptr<Loads::NodalLoadDataInfo> load_info
            (Loads::createNodalLoadDataInfo(NODAL_POINT_LOAD::num()).release());
          load_info->clear();
          load_info->setLoadCaseID(load_cases[i]);
          load_info->setLoadNameEnumID(NODAL_TEMPERATURE::num());
          load_info->setTime(time_vec.el(j));
          if (k>0)
            load_info->setDVID((*dv_vector)[k-1]->getID());
          
          unsigned int n_dofs = solution->size(),
            load_set_ID = load_database.getAvailableLoadSetID();
          unsigned int dv_ID = 0;
          if (k>0)
            dv_ID = (*dv_vector)[k-1]->getID();
          else
            dv_ID = FESystemNumbers::InvalidID;
          
          const_cast<TimeDependentLoadCase&>(load_case).addLoadSetID
            (time_vec.el(j), NODAL_TEMPERATURE::num(), load_set_ID, 1.0, dv_ID);

          // stream to which the data will be written
          std::stringstream load_data;
          
          load_data << "NODAL_LOAD_SET BEGIN " << std::endl
            << "NAME NODAL_TEMPERATURE_LOAD" << std::endl
            << "ID " << load_set_ID <<  std::endl
            << "KIND NODAL_TEMPERATURE" << std::endl
            << "N_DOFS 1" << std::endl
            << "N_LOADS " << n_dofs << std::endl;
          
          MeshBase::const_node_iterator           node  = from_mesh->nodes_begin();
          const MeshBase::const_node_iterator end_node  = from_mesh->nodes_end();
          
          load_ID = 1;
          for ( ; node != end_node; node++)
            {
            node_ID = from_mesh_data->getForeignIDFromNode(*node);
            dof_number = (*node)->dof_number(0,0,0);
            
            load_data << load_ID << "  " 
              << node_ID << "  " 
              << std::setprecision(16) << (*solution.get())(dof_number) 
              << std::endl;
            
            load_ID++;
            }
          
          load_data << "NODAL_LOAD_SET END" << std::endl;

          // now read and store the load data
          std::string tag;
          FESystemIOUtility::peekFromInput(load_data, tag);
          unsigned int enum_ID = LoadSetKindEnum::enumID(tag);
          std::auto_ptr<LoadSetBase> load_set(createLoadSet(enum_ID).release());
          load_set->readFromInputStream(load_data);
          load_database.addLoadSet(load_set.release());
          }
        }
      }
    }
  else
    {
    std::auto_ptr<FESystemDatabase::TimeIndependentDataInfo> data_info
    (new FESystemDatabase::TimeIndependentDataInfo);
    
    for (unsigned int i=0; i < load_cases.size(); i++) // iterate over the load cases
      {
      for (unsigned int k=0; k < (dv_vector->size() + 1); k++) // iterate over the load cases
        {
        // clear the data info and set the data for the solution vector 
        // to be read in
        data_info->clear();
        data_info->setDisciplineEnumID(from_discipline_enum_ID);
        data_info->setName("Solution");
        if (k>0) // write the dv value
          data_info->setDVID((*dv_vector)[k-1]->getID());
        data_info->setLoadCase(load_cases[i]);
        
        this->fesystem_controller.global_data_storage->fillVector(*data_info, *solution);
        
        // make sure that the number of dofs in the mesh are the same as the number of dofs 
        // in this solution vector
        Assert(solution->size() == from_dof_map->n_dofs(), ExcInternalError());
        
        // iterate over all the dofs, and write the load values to the output stream
        // create a load data info for the loads
        // get the load case from the database
        const TimeIndependentLoadCase& load_case = dynamic_cast<const TimeIndependentLoadCase&>
          (load_database.getLoadCaseFromID(load_cases[i]));
        
        std::auto_ptr<Loads::NodalLoadDataInfo> load_info
          (Loads::createNodalLoadDataInfo(NODAL_POINT_LOAD::num()).release());
        load_info->clear();
        load_info->setLoadCaseID(load_cases[i]);
        load_info->setLoadNameEnumID(NODAL_TEMPERATURE::num());
        if (k>0)
          load_info->setDVID((*dv_vector)[k-1]->getID());
        
        unsigned int n_dofs = solution->size(),
          load_set_ID = load_database.getAvailableLoadSetID();
        unsigned int dv_ID = 0;
        if (k>0)
          dv_ID = (*dv_vector)[k-1]->getID();
        else
          dv_ID = FESystemNumbers::InvalidID;
        
        
        const_cast<TimeIndependentLoadCase&>(load_case).addLoadSetID
          (NODAL_TEMPERATURE::num(), load_set_ID, 1.0, dv_ID);
        
        // stream to which the data will be written
        std::stringstream load_data;
        
        load_data << "NODAL_LOAD_SET BEGIN " << std::endl
          << "NAME NODAL_TEMPERATURE_LOAD" << std::endl
          << "ID " << load_set_ID <<  std::endl
          << "KIND NODAL_TEMPERATURE" << std::endl
          << "N_DOFS 1" << std::endl
          << "N_LOADS " << n_dofs << std::endl;
        
        MeshBase::const_node_iterator           node  = from_mesh->nodes_begin();
        const MeshBase::const_node_iterator end_node  = from_mesh->nodes_end();
        
        load_ID = 1;
        for ( ; node != end_node; node++)
          {
          node_ID = from_mesh_data->getForeignIDFromNode(*node);
          dof_number = (*node)->dof_number(0,0,0);
          
          load_data << load_ID << "  " 
            << node_ID << "  " 
            << std::setprecision(16) << (*solution.get())(dof_number) 
            << std::endl;
          
          load_ID++;
          }
        
        load_data << "NODAL_LOAD_SET END" << std::endl;
        
        // now read and store the load data
        std::string tag;
        FESystemIOUtility::peekFromInput(load_data, tag);
        unsigned int enum_ID = LoadSetKindEnum::enumID(tag);
        std::auto_ptr<LoadSetBase> load_set(createLoadSet(enum_ID).release());
        load_set->readFromInputStream(load_data);
        load_database.addLoadSet(load_set.release());
        }
      }
    }
}

