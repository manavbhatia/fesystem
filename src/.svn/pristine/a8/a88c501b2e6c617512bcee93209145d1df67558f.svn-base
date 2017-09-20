// $Id: GmshOutputProcessor.C,v 1.1.6.4 2008-08-21 00:55:20 manav Exp $

// C++ includes
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <sstream>


// FESystem includes
#include "OutputProcessor/GmshOutputProcessor.h"
#include "FESystem/FESystemController.h"
#include "FESystem/AnalysisCase.h"
#include "PostProcess/PostProcessQtyDatabase.h"
#include "PostProcess/ElemPostProcessQty.h"
#include "Mesh/FEMeshData.h"
#include "Mesh/MeshList.h"
#include "Database/GlobalDataStorage.h"
#include "Discipline/ThermalDisciplineInfo.h"
#include "Discipline/StructuralDisciplineInfo.h"
#include "Discipline/PistonTheory.h"
#include "OutputProcessor/OutputInfo.h"
#include "Solutions/SolutionBase.h"
#include "Utilities/AutoptrVector.h"


// libmesh includes
#include "geom/elem.h"
#include "geom/node.h"
#include "numerics/numeric_vector.h"


GmshOutputProcessor::GmshOutputProcessor(FESystem::FESystemController& fesys_controller,
                                         const OutputInfo& info):
OutputProcessor(fesys_controller, info, GMSH_OUTPUT_PROCESSOR::num()),
mesh(NULL),
mesh_data(NULL),
dummy_elem_set(NULL)
{
	
}




GmshOutputProcessor::~GmshOutputProcessor()
{
	
}



void GmshOutputProcessor::writeData()
{
  // iterate over the disciplines and write the mesh to 
  // each discipline's file
  
  std::vector<Discipline::DisciplineInfo*>& disciplines = 
  this->fesystem_controller.analysis_case->getAnalysisDisciplineInfos();
  
  unsigned int mesh_ID =0;
  
  std::vector<Discipline::DisciplineInfo*>::const_iterator it, end;
  it = disciplines.begin();
  end = disciplines.end();
  
  for (; it !=end ; it++)
    {
      mesh_ID = (*it)->getMeshID();
      
      // get the mesh and meshdata from mesh list
      this->mesh = 
      this->fesystem_controller.mesh_list->getMeshFromID(mesh_ID);
      this->mesh_data = 
      this->fesystem_controller.mesh_list->getMeshDataFromID(mesh_ID);
      this->dummy_elem_set = 
      this->fesystem_controller.mesh_list->getDummyElemSetFromID(mesh_ID);
      
      // get the name of the file to be written to
      std::string file_name = 
      Discipline::AnalysisDisciplineEnum::enumName((*it)->getDisciplineEnumID());
      
      std::string name = this->output_info.getOutputFileName();
      
      file_name += name;
      
      // open the file, and call the functions to write the data to this stream
      std::fstream output_file;
      output_file.open(file_name.c_str(), std::fstream::out);
      
      unsigned int solution_tensor_order = 0;
      // set the tensor order for the discipline
      switch ((*it)->getDisciplineEnumID())
      {
        case THERMAL_DISCIPLINE_ENUM_ID:
          solution_tensor_order = 0;
          break;
          
        case STRUCTURAL_DISCIPLINE_ENUM_ID:
          solution_tensor_order = 1;
          break;
          
        case PISTON_THEORY_ENUM_ID:
          continue;
          break;
          
        default:
          Assert(false, ExcInternalError());
      }
      
      this->writeNodalData(output_file, (*it)->getDisciplineEnumID(), 
                           solution_tensor_order);
      
      output_file.close();
    }
}



std::string
GmshOutputProcessor::getElemTypeName(ElemType type, 
                                     const unsigned int tensor_order)
{
  static std::string name1, name2;
  name1.clear(); name2.clear();
  
  switch(type)
  {
    case EDGE2:
      name1 = "L";
      break;
      
    case TRI3:
      name1 = "T";
      break;
      
    case QUAD4:
      name1 = "Q";
      break;
      
    case HEX8:
      name1 = "H";
      break;
      
    case PRISM6:
      name1 = "I";
      break;
      
    default:
      abort();
  }
  
  switch (tensor_order)
  {
    case 0:
      name2 = "S";
      break;
      
    case 1:
      name2 = "V";
      break;
      
    case 2:
      name2 = "T";
      break;
      
    default:
      abort();
  }
  
  name2 += name1;
  return name2;
}






void GmshOutputProcessor::writeNodalData(std::ostream& output,
                                         const unsigned int discipline,
                                         const unsigned int tensor_order)
{
  
  // prepare the elem type map
  this->prepareElemTypeMap();
  
  // write the file header
  output << "$PostFormat" << std::endl
  << 1.4 << "  " << 0 << "  " << sizeof(double) << std::endl
  << "$EndPostFormat" << std::endl;
  
  MeshBase::const_element_iterator elem_it = this->mesh->elements_begin();
  MeshBase::const_element_iterator elem_end = this->mesh->elements_end();
  
  std::set<Elem*>::const_iterator dummy_elem_end =
  this->dummy_elem_set->end();
  
  // iterate over each solution, and prepare a list of nodal and element centered 
  // variables for the same.
  const std::vector<Solution::SolutionBase*> sols = 
  this->fesystem_controller.analysis_case->getSolutions();
  
  std::vector<Solution::SolutionBase*>::const_iterator sol_it, sol_end;
  sol_it = sols.begin();
  sol_end = sols.end();
  
  std::vector<unsigned int>::const_iterator lc_it, lc_end;
  
  unsigned int n_dofs = pow(3,tensor_order);
  output.setf(std::ios::showpoint);
  
  // set the number of scalars, vectors and tensors for solutions
  unsigned int n_scalars = 0, n_vectors = 0, n_tensors = 0;
  
  if (tensor_order == 0)
    n_scalars = 1;
  else if (tensor_order == 1)
    n_vectors = 1;
  else if (tensor_order == 2)
    n_tensors = 1;
  else
    Assert(false, ExcInternalError());
  
  FESystemUtility::AutoPtrVector<NumericVector<double> > sol_vec;
  
  // iterate over all the solutions and write the nodal solution names for the 
  // load cases and DV sensitivities.
  for ( ; sol_it != sol_end; sol_it++)
    {
      // if the solution does not have this discipline, nothing to be written.
      if (! (*sol_it)->checkParticipatingDiscipline(discipline))
        continue;
      
      // write the time independent data
      // else, get the solution name and load cases from this solution
      FESystemUtility::AutoPtrVector<FESystemDatabase::DataInfoBase> 
      sol_names((*sol_it)->getTimeIndependentNodalSolutionDataInfo(discipline).release());
      
      sol_vec.clear();
      sol_vec.push_back(NumericVector<double>::build().release());
      
      for (unsigned int i=0; i < sol_names.size(); i++)
        {
          output << "$View " << std::endl
          // only one time step for a time independent solution
          << sol_names[i]->getDisplayName() <<  " " << 1 << std::endl;
          
          
          this->writeViewHeader(output, n_scalars, n_vectors, n_tensors);
          
          // static solution, hence, only at the first time step
          output << 0.0 << std::endl; 
          
          this->fesystem_controller.global_data_storage->fillVector(*sol_names[i],
                                                                    *(sol_vec[0]));
          
          this->writeSolutionVector(output, sol_vec.get(), n_dofs);
          
          output << "$EndView" << std::endl;
        }
      
      // delete the solution vector that was used for the static solution
      sol_vec.clear();
      
      // write the time dependent data
      // get the solution name and load cases from this solution
      FESystemUtility::AutoPtrVector<FESystemDatabase::TransientDataInfoSet >
      time_dep_sol_sets
      ((*sol_it)->getTimeDependentNodalSolutionDataInfo(discipline).release());
      
      unsigned int n_data_info = 0;
      
      for (unsigned int i=0 ; i < time_dep_sol_sets.size(); i++)
        {
          n_data_info = time_dep_sol_sets[i]->getNDataInfo();
          if (n_data_info == 0)
            continue;
          
          FESystemUtility::AutoPtrVector<FESystemDatabase::TimeDependentDataInfo> 
          time_data_info(time_dep_sol_sets[i]->getDataInfoForAllTimeSteps().release());
          
          // add the solutions to the vec
          if (sol_vec.size() != n_data_info)
            {
              sol_vec.clear();
              sol_vec.resize(n_data_info);
              for (unsigned int j=0; j < n_data_info; j++)
                {
                  sol_vec.reset(j, NumericVector<double>::build().release());
                  this->fesystem_controller.global_data_storage->fillVector
                  (*(time_data_info[j]), *(sol_vec[j]));
                }
            }
          
          // write the name and the number of time sets
          output << "$View " << std::endl
          // the display name will be taken from the first solution
          << time_dep_sol_sets[i]->getDisplayName() <<  " "
          << n_data_info << std::endl;
          
          this->writeViewHeader(output, n_scalars, n_vectors, n_tensors);
          
          for (unsigned int j=0; j < n_data_info; j++)
            output << time_data_info[j]->getTimeValue() << "  ";
          
          this->writeSolutionVector(output, sol_vec.get(), n_dofs);
          
          output << "$EndView" << std::endl;
        }
    }
  output << std::endl;
}




void GmshOutputProcessor::prepareElemTypeMap()
{
  // first of all, clear the map
  this->elem_vec_map.clear();
  
  // next, iterate over all the elems, and based on the elem type, 
  // insert them in the vector. Finally, insert the vectors in the 
  // map
  // get the iterator for the elems
  MeshBase::const_element_iterator elem_it = this->mesh->elements_begin();
  MeshBase::const_element_iterator elem_end = this->mesh->elements_end();
  
  std::set<Elem*>::const_iterator dummy_elem_end =
  this->dummy_elem_set->end();
  
  std::vector<Elem*> *elem_vec=NULL;
	
  for (; elem_it != elem_end; elem_it++)
    {
      if (this->dummy_elem_set->find(*elem_it) != dummy_elem_end)
        continue;
      
      
      // insert this elem type in the map if it does not already exist
      if (this->elem_vec_map.count((**elem_it).type()) == 0)
        {
          elem_vec = 
          &(this->elem_vec_map.insert(std::map<ElemType, std::vector<Elem*> >::value_type
                                      ((**elem_it).type(), std::vector<Elem*>())).first->second);
        }
      else
        elem_vec = &(this->elem_vec_map.find((**elem_it).type())->second);
      
      elem_vec->push_back(*elem_it);
    }
  
}



void GmshOutputProcessor::writeViewHeader(std::ostream& output, 
                                          const unsigned int n_scalars,
                                          const unsigned int n_vectors,
                                          const unsigned int n_tensors)
{
  // write the number of scalar, vector, tensors 
  unsigned int n_entity;
  
  // first for the points
  output << 0 << " " << 0 <<  " " <<  0 << std::endl;
  // first order entities
  
  // lines
  n_entity = 
  (this->elem_vec_map.count(EDGE2)==0)?0:this->elem_vec_map.find(EDGE2)->second.size();
  output << n_scalars * n_entity << "  "
  << n_vectors * n_entity << "  " 
  << n_tensors * n_entity << "  " << std::endl;
  
  // triangles
  n_entity = 
  (this->elem_vec_map.count(TRI3)==0)?0:this->elem_vec_map.find(TRI3)->second.size();
  output << n_scalars * n_entity << "  "
  << n_vectors * n_entity << "  " 
  << n_tensors * n_entity << "  " << std::endl;
  
  // quadrangles
  n_entity = 
  (this->elem_vec_map.count(QUAD4)==0)?0:this->elem_vec_map.find(QUAD4)->second.size();
  output << n_scalars * n_entity << "  "
  << n_vectors * n_entity << "  " 
  << n_tensors * n_entity << "  " << std::endl;
  
  // tetrahedra
  n_entity = 
  (this->elem_vec_map.count(TET4)==0)?0:this->elem_vec_map.find(TET4)->second.size();
  output << n_scalars * n_entity << "  "
  << n_vectors * n_entity << "  " 
  << n_tensors * n_entity << "  " << std::endl;
  
  // hexahedra
  n_entity = 
  (this->elem_vec_map.count(HEX8)==0)?0:this->elem_vec_map.find(HEX8)->second.size();
  output << n_scalars * n_entity << "  "
  << n_vectors * n_entity << "  " 
  << n_tensors * n_entity << "  " << std::endl;
  
  // prisms
  n_entity = 
  (this->elem_vec_map.count(PRISM6)==0)?0:this->elem_vec_map.find(PRISM6)->second.size();
  output << n_scalars * n_entity << "  "
  << n_vectors * n_entity << "  " 
  << n_tensors * n_entity << "  " << std::endl;
  
  // pyramids
  n_entity = 
  (this->elem_vec_map.count(PYRAMID5)==0)?0:this->elem_vec_map.find(PYRAMID5)->second.size();
  output << n_scalars * n_entity << "  "
  << n_vectors * n_entity << "  " 
  << n_tensors * n_entity << "  " << std::endl;
  
  
  // second order entities
  
  // lines
  n_entity = 
  (this->elem_vec_map.count(EDGE3)==0)?0:this->elem_vec_map.find(EDGE3)->second.size();
  // the edge 4 will be written as a 2nd order element
  n_entity += 
  (this->elem_vec_map.count(EDGE4)==0)?0:this->elem_vec_map.find(EDGE4)->second.size();
  output << n_scalars * n_entity << "  "
  << n_vectors * n_entity << "  " 
  << n_tensors * n_entity << "  " << std::endl;
  
  // triangles
  n_entity = 
  (this->elem_vec_map.count(TRI6)==0)?0:this->elem_vec_map.find(TRI6)->second.size();
  output << n_scalars * n_entity << "  "
  << n_vectors * n_entity << "  " 
  << n_tensors * n_entity << "  " << std::endl;
  
  // quadrangles
  n_entity = 
  (this->elem_vec_map.count(QUAD8)==0)?0:this->elem_vec_map.find(QUAD8)->second.size();
  n_entity += 
  (this->elem_vec_map.count(QUAD9)==0)?0:this->elem_vec_map.find(QUAD9)->second.size();
  output << n_scalars * n_entity << "  "
  << n_vectors * n_entity << "  " 
  << n_tensors * n_entity << "  " << std::endl;
  
  // tetrahedra
  n_entity = 
  (this->elem_vec_map.count(TET10)==0)?0:this->elem_vec_map.find(TET10)->second.size();
  output << n_scalars * n_entity << "  "
  << n_vectors * n_entity << "  " 
  << n_tensors * n_entity << "  " << std::endl;
  
  // hexahedra
  n_entity = 
  (this->elem_vec_map.count(HEX20)==0)?0:this->elem_vec_map.find(HEX20)->second.size();
  n_entity += 
  (this->elem_vec_map.count(HEX27)==0)?0:this->elem_vec_map.find(HEX27)->second.size();
  output << n_scalars * n_entity << "  "
  << n_vectors * n_entity << "  " 
  << n_tensors * n_entity << "  " << std::endl;
  
  // prisms
  n_entity = 
  (this->elem_vec_map.count(PRISM15)==0)?0:this->elem_vec_map.find(PRISM15)->second.size();
  n_entity += 
  (this->elem_vec_map.count(PRISM18)==0)?0:this->elem_vec_map.find(PRISM18)->second.size();
  output << n_scalars * n_entity << "  "
  << n_vectors * n_entity << "  " 
  << n_tensors * n_entity << "  " << std::endl;
  
  // pyramids
  n_entity = 0;
  output << n_scalars * n_entity << "  "
  << n_vectors * n_entity << "  " 
  << n_tensors * n_entity << "  " << std::endl;
  
  // finally, texts
  output << "0 0 0 0 " << std::endl;
  
}



void GmshOutputProcessor::writeSolutionVector(std::ostream& output,
                                              const std::vector<NumericVector<double>*>* sol_vec,
                                              const unsigned int n_dofs)
{
  
  // a vector is created to decide on how to iterate over the elements
  std::vector<ElemType> elem_type_vec;
  elem_type_vec.push_back(EDGE2);
  elem_type_vec.push_back(TRI3);
  elem_type_vec.push_back(QUAD4);
  elem_type_vec.push_back(TET4);
  elem_type_vec.push_back(HEX8);
  elem_type_vec.push_back(PRISM6);
  elem_type_vec.push_back(PYRAMID5);
  // the second order elements are being left out for now, since the code does not have
  // support for that
  
  unsigned int n_nodes=0, dof = 0;
  std::vector<ElemType>::const_iterator elem_type_it, elem_type_end;
  elem_type_it = elem_type_vec.begin();
  elem_type_end = elem_type_vec.end();
  
  for (; elem_type_it != elem_type_end; elem_type_it++)
    {
      // get the vector of elemets from the map
      if (this->elem_vec_map.count(*elem_type_it) == 0)
        continue;
      
      const std::vector<Elem*>& elem_vec = this->elem_vec_map.find(*elem_type_it)->second;
      std::vector<Elem*>::const_iterator elem_it, elem_end;
      elem_it = elem_vec.begin();
      elem_end = elem_vec.end();
      
      // iterate over the elements and write them
      for (; elem_it != elem_end; elem_it++)
        {
          n_nodes = (**elem_it).n_nodes();
          
          for (unsigned int j=0; j < 3; j++)
            {
              for (unsigned int i=0; i < n_nodes; i++)
                output << (**elem_it).point(i)(j) << " ";
              output << std::endl;
            }
          output << std::endl;
          
          // iterate over all the solutions, and write them
          for (unsigned int k=0; k < sol_vec->size(); k++)
            {
              for (unsigned int i=0; i < n_nodes; i++)
                {
                  for (unsigned int j=0; j < n_dofs; j++)
                    {
                      dof = (**elem_it).get_node(i)->dof_number(0,j,0);
                      output <<  std::setprecision(8) << (*(*sol_vec)[k])(dof) << " ";
                    }
                  output << std::endl;
                }
              output << std::endl;
            }
          output << std::endl;
        }
    }
}

