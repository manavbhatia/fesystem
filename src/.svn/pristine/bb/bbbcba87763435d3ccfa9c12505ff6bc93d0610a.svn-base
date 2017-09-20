// $Id: InterpolationDriver.C,v 1.12 2006-10-23 23:42:35 manav Exp $

// C++ includes
#include <sstream>
#include <memory>

// FESystem includes
#include "Interpolation/InterpolationDriver.h"
#include "Interpolation/InterpolationCase.h"
#include "FESystem/FESystemController.h"
#include "FESystem/AnalysisCase.h"
#include "Mesh/MeshList.h"
#include "Mesh/FEMeshData.h"
#include "Database/GlobalDataStorage.h"
#include "Loads/LoadDatabase.h"
#include "Interpolation/FEInterpolation.h"
#include "DesignData/DesignDatabase.h"
#include "Discipline/StructuralAnalysis.h"


// libmesh includes
#include "base/dof_map.h"
#include "mesh/mesh.h"
#include "base/node.h"
#include "numerics/numeric_vector.h"



InterpolationDriver::InterpolationDriver(FESystem::FESystemController& controller,
                                         InterpolationCase& interp_case):
fesystem_controller(controller),
interpolation_case(interp_case),
interpolation_base(NULL),
mesh(NULL),
mesh_data(NULL),
dof_map(NULL)
{
  this->init();
}





InterpolationDriver::~InterpolationDriver()
{
  
}





void InterpolationDriver::init()
{
  // the analysis here will proceede as follows
  // check the type of interpolation requested
  switch (this->interpolation_case.type())
    {
    case InterpolationCase::DIRECT:
      {
        
      }
      break;
      
    case InterpolationCase::FE:
      {
        this->interpolation_base.reset(new FEInterpolation(*this));
      } 
      break;
      
    case InterpolationCase::LEAST_SQUARE:
      {
        //this->interpolation_base.reset
        //  (new LSInterpolation);
      }
      break;
      
    case InterpolationCase::SPLINE:
      {
        //	this->interpolation_base.reset
        //(new SplineInterpolation);
      }
      break;
      
    default:
      abort();
      break;
    }
  
  // init the mesh, mesh data and the dof_map for the 
  // related structures, with an exception of the direct interpolation
  if (this->interpolation_case.type() == InterpolationCase::DIRECT)
    this->mesh_ID = 
      this->interpolation_case.fromMeshID();
  else 
    this->mesh_ID = 
      this->interpolation_case.toMeshID();
  
  // now init the data structures
  this->mesh = 
    this->fesystem_controller.mesh_list->getMeshFromID( this->mesh_ID);
  this->mesh_data = 
    this->fesystem_controller.mesh_list->getMeshDataFromID( this->mesh_ID);
  this->dof_map = 
    this->fesystem_controller.mesh_list->getDofMapFromID( this->mesh_ID);
  
  
}



void InterpolationDriver::interpolateAndCreateLoads(const std::vector<unsigned int>& load_cases)
{
  // if this is not direct transfer of loads, 
  // interpolate, else just create loads
  if (this->interpolation_case.type() != InterpolationCase::DIRECT)
    {
    this->interpolate(load_cases);
    }
  
  
  // create loads after interpolation is done
  this->createLoads(load_cases);
}




void InterpolationDriver::interpolate(const std::vector<unsigned int>& load_cases)
{
  std::auto_ptr<NumericVector<double> > 
  solution(NumericVector<double>::build().release()),
  interpolated_solution(NumericVector<double>::build().release());
  
  std::string sol_name = "Thermal_Solution",
  interpolated_sol_name = "InterpolatedThermal_Solution";
  
  GlobalDataStorage& global_data_storage = 
    *(this->fesystem_controller.global_data_storage.get());
  
  DofMap& to_mesh_dof_map = 
    *(this->fesystem_controller.mesh_list->getDofMapFromID
      (this->interpolation_case.toMeshID()));
  
  
//  // iterate over the load cases, and the design variables
//  // get the load case vector from the analysis load case
//  const std::vector<unsigned int>& load_cases = 
//    this->fesystem_controller.analysis_case->getLoadCaseIDs();
	
  std::vector<unsigned int>::const_iterator load_case_it, load_case_begin,
    load_case_end;
  
  load_case_begin = load_cases.begin();
  load_case_end = load_cases.end();	
  load_case_it = load_case_begin;
  
  // iterate over each load case, and get the interpolated vector from the base
  for (; load_case_it != load_case_end; load_case_it++)
    {
    // get the solution in the vector, and ask the interpolation_base 
    // to return the interpolated solution
    global_data_storage.fillVector(*load_case_it, 
                                   sol_name,
                                   *(solution.get()));
    if (!interpolated_solution->initialized() ||
        interpolated_solution->size() != to_mesh_dof_map.n_dofs())
      interpolated_solution->init(to_mesh_dof_map.n_dofs());
    
    this->interpolation_base->getInterpolatedValues((*solution.get()),
                                                    (*interpolated_solution.get()));
    
    // save the interpolated value
    this->fesystem_controller.global_data_storage->storeVector
      (*load_case_it, 
       interpolated_sol_name, 
       (*interpolated_solution.get()));
    }
  
  
  // now iterate over the design variables
  load_case_it = load_case_begin;
  
  std::auto_ptr<std::vector<DesignData::DesignParameter*> > dv_vector =
    this->fesystem_controller.design_database->getParameters();
  
  std::vector<DesignData::DesignParameter*>::const_iterator dv_it, dv_end;
  dv_end = dv_vector->end();
  
  std::string sol_sens_name, interpolated_sol_sens_name;
  
  // iterate over each load case, and ask solver to solve for it
  for (; load_case_it != load_case_end; load_case_it++)
    {
    dv_it = dv_vector->begin();
		
    for (; dv_it != dv_end; dv_it++)
      {
      // create the solution name
      std::ostringstream dv_ID;
      dv_ID << (*dv_it)->getID();
      sol_sens_name.clear();
      sol_sens_name = "d";
      sol_sens_name += sol_name;
      sol_sens_name += "_dDV";
      sol_sens_name += dv_ID.str();
      
      interpolated_sol_sens_name.clear();
      interpolated_sol_sens_name = "d";
      interpolated_sol_sens_name += interpolated_sol_name;
      interpolated_sol_sens_name += "_dDV";
      interpolated_sol_sens_name += dv_ID.str();
      
      
      // get the solution in the vector, and ask the interpolation_base 
      // to return the interpolated solution
      global_data_storage.fillVector(*load_case_it, 
                                     sol_sens_name,
                                     *(solution.get()));
      if (interpolated_solution->size() != to_mesh_dof_map.n_dofs())
        interpolated_solution->init(to_mesh_dof_map.n_dofs());
      
      this->interpolation_base->getInterpolatedValues((*solution.get()),
                                                      (*interpolated_solution.get()));
      
      // save the interpolated value
      this->fesystem_controller.global_data_storage->storeVector
        (*load_case_it, 
         interpolated_sol_sens_name, 
         (*interpolated_solution.get()));
      
      }
    }
}



void InterpolationDriver::createLoads(const std::vector<unsigned int>& load_cases)
{
  unsigned int load_ID = 0, load_set_ID = 0, n_nodes = 0,
  node_ID = 0, dof_number = 0, mesh_ID = 0;
  std::string interpolated_sol_name, interpolated_sol_sens_name,
  load_sens_name;
  
  
  // reference to load database and global database
  LoadDatabase& load_database = 
    *(this->fesystem_controller.load_database.get());
  GlobalDataStorage& global_data_storage = 
    *(this->fesystem_controller.global_data_storage.get());
  MeshDS::MeshList& mesh_list = 
    *(this->fesystem_controller.mesh_list.get());
  
  std::auto_ptr<NumericVector<double> > 
    solution(NumericVector<double>::build().release());
  
  // get the mesh and mesh data for the mesh to which the 
  // data has been interpolated
  if (this->interpolation_case.type() == InterpolationCase::DIRECT)
    {
    interpolated_sol_name = "Thermal_Solution";
    mesh_ID = this->interpolation_case.fromMeshID();
    }
  else 
    {
    interpolated_sol_name = "InterpolatedThermal_Solution";
    mesh_ID = this->interpolation_case.toMeshID();
    }
  
  
  MeshDS::FEMesh* mesh = mesh_list.getMeshFromID(mesh_ID);
  MeshDS::FEMeshData* mesh_data = mesh_list.getMeshDataFromID(mesh_ID);
  
  
//  // iterate over the load cases, and the design variables
//  // get the load case vector from the analysis load case
//  const std::vector<unsigned int>& load_cases = 
//    this->fesystem_controller.analysis_case->getLoadCaseIDs();	
  
  
  std::vector<unsigned int>::const_iterator load_case_it, load_case_begin,
    load_case_end;
  
  load_case_begin = load_cases.begin();
  load_case_end = load_cases.end();	
  load_case_it = load_case_begin;
  
  
  // iterate over each load case, and get the interpolated vector from the base
  for (; load_case_it != load_case_end; load_case_it++)
    {
    // read the load vector
    global_data_storage.fillVector(*load_case_it, 
                                   interpolated_sol_name,
                                   *(solution.get()));
    n_nodes = solution->size();
    
    // get the load case ID for this load
    // for now, this is being done for steady state
    load_set_ID = load_database.getLoadSetID(*load_case_it,
                                             0.0,
                                             NODAL_TEMPERATURE::num());
    
    std::stringstream load_data;
    
    // load set will begin with a tag and the number of nodes in it
    load_data << "BEGIN_LOAD_SET  " << n_nodes << std::endl;
    load_data << "LOAD_SET_NAME  " << NODAL_TEMPERATURE::name() << std::endl;
    load_data << "LOAD_SET_ID  " << load_set_ID << std::endl;
    load_data << "LOAD_SET_KIND  POINT_LOAD " << std::endl;
    load_data << "LOAD_SET_N_DOFS  1 " << std::endl;
    
    // iterate over each node in the mesh to which the 
    // data is being interpolated
    load_ID = 1;
    // iterate over all nodes
    MeshBase::const_node_iterator           node  = mesh->nodes_begin();
    const MeshBase::const_node_iterator end_node  = mesh->nodes_end();
    
    for ( ; node != end_node; ++node)
      {
      node_ID = mesh_data->getForeignIDFromNode(*node);
      dof_number = (*node)->dof_number(0,0,0);
      
      load_data << load_ID << "  " 
		    << node_ID << "  " 
		    << std::setprecision(16) << (*solution.get())(dof_number) 
		    << std::endl;
      
      load_ID++;
      }
    
    load_data << "END_LOAD_SET " << std::endl;
    
    // next, create a LoadSet, read it from the iostream, 
    // and add it to the load_database
    std::auto_ptr<LoadSet> load_set(new LoadSet);
    load_set->readFromInputStream(load_data);
    load_ID = load_set->getLoadSetID();
    load_database.addLoadSet(load_ID,
                             load_set.release());
    }
  
  
  
  
  // now iterate over the design variables
  load_case_it = load_case_begin;
  
  std::auto_ptr<std::vector<DesignData::DesignParameter*> > dv_vector =
    this->fesystem_controller.design_database->getParameters();
  
  std::vector<DesignData::DesignParameter*>::const_iterator dv_it, dv_end;
  dv_end = dv_vector->end();
	
  // iterate over each load case, and ask solver to solve for it
  for (; load_case_it != load_case_end; load_case_it++)
    {
    dv_it = dv_vector->begin();
		
    for (; dv_it != dv_end; dv_it++)
      {
      std::ostringstream dv_ID;
      dv_ID << (*dv_it)->getID();
      
      // read the solution
      interpolated_sol_sens_name.clear();
      interpolated_sol_sens_name = "Thermal_dSolution";
      interpolated_sol_sens_name += "_dDV_";
      interpolated_sol_sens_name += dv_ID.str();
      
      
      global_data_storage.fillVector(*load_case_it, 
                                     interpolated_sol_sens_name,
                                     *(solution.get()));
      n_nodes = solution->size();
      
      // get the load case ID for this load
      load_set_ID = load_database.getLoadSetIDForLoadSensitivity
        (*load_case_it, (*dv_it)->getID(), 0.0, NODAL_TEMPERATURE::num());
      
      std::stringstream load_data;
      
      load_sens_name.clear();
      load_sens_name = NODAL_TEMPERATURE::name();
      load_sens_name += "_Sens";
      
      // load set will begin with a tag and the number of nodes in it
      load_data << "BEGIN_LOAD_SET  " << n_nodes << std::endl;
      load_data << "LOAD_SET_NAME  " << load_sens_name << std::endl;
      load_data << "LOAD_SET_ID  " << load_set_ID << std::endl;
      load_data << "LOAD_SET_KIND  POINT_LOAD " << std::endl;
      load_data << "LOAD_SET_N_DOFS  1 " << std::endl;
      
      // iterate over each node in the mesh to which the 
      // data is being interpolated
      load_ID = 1;
      // iterate over all nodes
      MeshBase::const_node_iterator           node  = mesh->nodes_begin();
      const MeshBase::const_node_iterator end_node  = mesh->nodes_end();
      
      for ( ; node != end_node; ++node)
        {
	      node_ID = mesh_data->getForeignIDFromNode(*node);
	      dof_number = (*node)->dof_number(0,0,0);
        
	      load_data << load_ID << "  " 
          << node_ID << "  " 
          << std::setprecision(16) << (*solution.get())(dof_number) 
          << std::endl;
        
	      load_ID++;
        }
      
      load_data << "END_LOAD_SET " << std::endl;
      
      // next, create a LoadSet, read it from the iostream, 
      // and add it to the load_database
      std::auto_ptr<LoadSet> load_set(new LoadSet);
      load_set->readFromInputStream(load_data);
      load_ID = load_set->getLoadSetID();
      load_database.addLoadSet(load_ID,
                               load_set.release());
      
      }
    }
  
}
