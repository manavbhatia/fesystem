// $Id: InterpolationBase.C,v 1.5.4.2 2007-06-13 14:57:20 manav Exp $

// C++ includes

// FESystem includes
#include "Interpolation/InterpolationBase.h"
#include "Interpolation/InterpolationCase.h"
#include "Utilities/ElemSetList.h"
#include "FESystem/FESystemController.h"
#include "Mesh/MeshList.h"
#include "Mesh/FEMeshData.h"

// libMesh includes
#include "geom/point.h"
#include "geom/elem.h"
#include "numerics/dense_vector.h"

InterpolationBase::InterpolationBase(FESystem::FESystemController& controller, 
                                     const InterpolationCase& interp_case):
fesystem_controller(controller),
interpolation_case(interp_case)
//mesh_ID(FESystemNumbers::InvalidID),
//mesh(NULL),
//mesh_data(NULL),
//dof_map(NULL),
//sol_to_interpolate(NULL)
{
//  this->init(); 
}


InterpolationBase::~InterpolationBase()
{
  
}


//void
//InterpolationBase::init()
//{
//  // the analysis here will proceede as follows
//  // check the type of interpolation requested
//  switch (this->interpolation_case.type())
//    {
//    case InterpolationCase::DIRECT:
//      {
//        
//      }
//      break;
//      
//    case InterpolationCase::FE:
//      {
//        this->interpolation_base.reset(new FEInterpolation(*this));
//      } 
//      break;
//      
//    case InterpolationCase::LEAST_SQUARE:
//      {
//        //this->interpolation_base.reset
//        //  (new LSInterpolation);
//      }
//      break;
//      
//    case InterpolationCase::SPLINE:
//      {
//        //	this->interpolation_base.reset
//        //(new SplineInterpolation);
//      }
//      break;
//      
//    default:
//      abort();
//      break;
//    }
//  
//  // init the mesh, mesh data and the dof_map for the 
//  // related structures, with an exception of the direct interpolation
//  if (this->interpolation_case.type() == InterpolationCase::DIRECT)
//    this->mesh_ID = 
//      this->interpolation_case.fromMeshID();
//  else 
//    this->mesh_ID = 
//      this->interpolation_case.toMeshID();
//  
//  // now init the data structures
//  this->mesh = 
//    this->fesystem_controller.mesh_list->getMeshFromID( this->mesh_ID);
//  this->mesh_data = 
//    this->fesystem_controller.mesh_list->getMeshDataFromID( this->mesh_ID);
//  this->dof_map = 
//    this->fesystem_controller.mesh_list->getDofMapFromID( this->mesh_ID);
//  
//  
//}




//void
//InterpolationBase::interpolateAndCreateLoads(const std::vector<unsigned int>& load_cases)
//{
//  // if this is not direct transfer of loads, 
//  // interpolate, else just create loads
//  if (this->interpolation_case.type() != InterpolationCase::DIRECT)
//    {
//    this->interpolate(load_cases);
//    }
//  
//  
//  // create loads after interpolation is done
//  this->createLoads(load_cases);
//}




//void InterpolationBase::interpolate(const std::vector<unsigned int>& load_cases)
//{
//  std::auto_ptr<NumericVector<double> > 
//  solution(NumericVector<double>::build().release()),
//  interpolated_solution(NumericVector<double>::build().release());
//  
//  std::string sol_name = "Thermal_Solution",
//  interpolated_sol_name = "InterpolatedThermal_Solution";
//  
//  GlobalDataStorage& global_data_storage = 
//    *(this->fesystem_controller.global_data_storage.get());
//  
//  DofMap& to_mesh_dof_map = 
//    *(this->fesystem_controller.mesh_list->getDofMapFromID
//      (this->interpolation_case.toMeshID()));
//  
//  
////  // iterate over the load cases, and the design variables
////  // get the load case vector from the analysis load case
////  const std::vector<unsigned int>& load_cases = 
////    this->fesystem_controller.analysis_case->getLoadCaseIDs();
//	
//  std::vector<unsigned int>::const_iterator load_case_it, load_case_begin,
//    load_case_end;
//  
//  load_case_begin = load_cases.begin();
//  load_case_end = load_cases.end();	
//  load_case_it = load_case_begin;
//  
//  // iterate over each load case, and get the interpolated vector from the base
//  for (; load_case_it != load_case_end; load_case_it++)
//    {
//    // get the solution in the vector, and ask the interpolation_base 
//    // to return the interpolated solution
//    global_data_storage.fillVector(*load_case_it, 
//                                   sol_name,
//                                   *(solution.get()));
//    if (!interpolated_solution->initialized() ||
//        interpolated_solution->size() != to_mesh_dof_map.n_dofs())
//      interpolated_solution->init(to_mesh_dof_map.n_dofs());
//    
//    this->interpolation_base->getInterpolatedValues((*solution.get()),
//                                                    (*interpolated_solution.get()));
//    
//    // save the interpolated value
//    this->fesystem_controller.global_data_storage->storeVector
//      (*load_case_it, 
//       interpolated_sol_name, 
//       (*interpolated_solution.get()));
//    }
//  
//  
//  // now iterate over the design variables
//  load_case_it = load_case_begin;
//  
//  std::auto_ptr<std::vector<DesignData::DesignParameter*> > dv_vector =
//    this->fesystem_controller.design_database->getParameters();
//  
//  std::vector<DesignData::DesignParameter*>::const_iterator dv_it, dv_end;
//  dv_end = dv_vector->end();
//  
//  std::string sol_sens_name, interpolated_sol_sens_name;
//  
//  // iterate over each load case, and ask solver to solve for it
//  for (; load_case_it != load_case_end; load_case_it++)
//    {
//    dv_it = dv_vector->begin();
//		
//    for (; dv_it != dv_end; dv_it++)
//      {
//      // create the solution name
//      std::ostringstream dv_ID;
//      dv_ID << (*dv_it)->getID();
//      sol_sens_name.clear();
//      sol_sens_name = "d";
//      sol_sens_name += sol_name;
//      sol_sens_name += "_dDV";
//      sol_sens_name += dv_ID.str();
//      
//      interpolated_sol_sens_name.clear();
//      interpolated_sol_sens_name = "d";
//      interpolated_sol_sens_name += interpolated_sol_name;
//      interpolated_sol_sens_name += "_dDV";
//      interpolated_sol_sens_name += dv_ID.str();
//      
//      
//      // get the solution in the vector, and ask the interpolation_base 
//      // to return the interpolated solution
//      global_data_storage.fillVector(*load_case_it, 
//                                     sol_sens_name,
//                                     *(solution.get()));
//      if (interpolated_solution->size() != to_mesh_dof_map.n_dofs())
//        interpolated_solution->init(to_mesh_dof_map.n_dofs());
//      
//      this->interpolation_base->getInterpolatedValues((*solution.get()),
//                                                      (*interpolated_solution.get()));
//      
//      // save the interpolated value
//      this->fesystem_controller.global_data_storage->storeVector
//        (*load_case_it, 
//         interpolated_sol_sens_name, 
//         (*interpolated_solution.get()));
//      
//      }
//    }
//}
//


//void InterpolationBase::createLoads(const std::vector<unsigned int>& load_cases)
//{
//}



//Elem* InterpolationBase::getElemContainingPoint(Point& point, 
//                                                unsigned int elem_ID,
//                                                unsigned int mesh_ID)
//{
//  // get the elem set pairs for the interpolation case
//  InterpolationCase& interpolation_case = 
//  this->interpolation_driver.getInterpolationCase();
//  
//  ElemSetList& set_list = 
//    *(this->interpolation_driver.getFESystemController().elem_set_list.get());
//  
//  unsigned int from_mesh_ID = 
//    interpolation_case.fromMeshID();
//  
//  MeshDS::FEMeshData* from_mesh_data = 
//    this->interpolation_driver.getFESystemController().mesh_list->getMeshDataFromID(from_mesh_ID);
//  
//  // iterate over each pair, get the set for the second ID in the pair, 
//  // check if the set contains this elem ID
//  std::vector<std::pair<unsigned int, unsigned int> >& elem_set_pair = 
//    interpolation_case.getElemSetPairs();
//  
//  std::vector<std::pair<unsigned int, unsigned int> >::const_iterator it, end;
//  it = elem_set_pair.begin();
//  end = elem_set_pair.end();
//  
//  unsigned int to_set_ID=0, from_set_ID=0;
//  
//  for (; it != end; it++)
//    {
//    to_set_ID = it->second;
//    
//    const ElemSet& elem_set = 
//      set_list.getElemSetFromID(to_set_ID);
//    
//    // make sure that the element and mesh belong to this set
//    assert (elem_set.getMeshID() == mesh_ID);
//    assert (elem_set.getMeshID() == 
//            this->interpolation_driver.getInterpolationMeshID());
//    
//    // if this set does not contain the element, check the next pair
//    if (!elem_set.containsElem(elem_ID, mesh_ID))
//      continue;
//    
//    // now that the elem is in this set, get the from elem set
//    from_set_ID = it->first;
//    
//    const ElemSet& from_set =
//      set_list.getElemSetFromID(from_set_ID);
//    
//    // make sure that the mesh for this set is the mesh ID from which data
//    // is being interpolated
//    assert (from_set.getMeshID() == from_mesh_ID);
//    
//    // iterate over the elems in this set, get the elem from mesh data and
//    // check if the point lies in the element
//    Elem* from_mesh_elem = NULL;
//    const std::set<unsigned int>& elem_ID_set = 
//      from_set.getElemIDs();
//    std::set<unsigned int>::const_iterator elem_it, elem_end;
//    elem_it = elem_ID_set.begin();
//    elem_end = elem_ID_set.end();
//    
//    for ( ; elem_it != elem_end; elem_it++)
//      {
//      from_mesh_elem = 
//	    const_cast<Elem*>(from_mesh_data->getElemFromForeignID(*elem_it));
//      
//      if (from_mesh_elem->contains_point(point))
//        return from_mesh_elem;
//      }
//    }
//  
//  // if the execution reaches here, then that implies that all sets and elems 
//  // have been checked, and none of the elements contains this point. This is an error,
//  // hence, abort
//  
//  abort();
//  return NULL;
//}
//

//void InterpolationBase::getDofValuesForSourceElem(Elem* elem,
//                                                  DenseVector<double>& vector)
//{
//  unsigned int n_nodes = elem->n_nodes();
//  
//  if (vector.size() != n_nodes)
//    vector.resize(n_nodes);
//  
//  vector.zero();
//  unsigned int dof_id = 0;
//  for (unsigned int i=0; i < n_nodes; i++)
//    {
//    dof_id = elem->get_node(i)->dof_number(0,0,0);
//    vector(i) = (*this->sol_to_interpolate)(dof_id);
//    }
//}
