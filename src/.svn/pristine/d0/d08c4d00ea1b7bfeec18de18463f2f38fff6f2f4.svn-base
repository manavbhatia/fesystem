// $Id: TecPlotOutputProcessor.C,v 1.15.6.2 2008-08-21 00:55:21 manav Exp $

// C++ includes
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <sstream>


// FESystem includes
#include "OutputProcessor/TecPlotOutputProcessor.h"
#include "FESystem/FESystemController.h"
#include "FESystem/AnalysisCase.h"
#include "PostProcess/PostProcessQtyDatabase.h"
#include "PostProcess/ElemPostProcessQty.h"
#include "Mesh/FEMeshData.h"
#include "Mesh/MeshList.h"
#include "Database/GlobalDataStorage.h"
#include "Discipline/ThermalDisciplineInfo.h"
#include "Discipline/StructuralDisciplineInfo.h"
#include "OutputProcessor/OutputInfo.h"
#include "Solutions/SolutionBase.h"


// libmesh includes
#include "geom/elem.h"
#include "geom/node.h"
#include "numerics/numeric_vector.h"


TecPlotOutputProcessor::TecPlotOutputProcessor(FESystem::FESystemController& fesys_controller,
                                               const OutputInfo& info):
OutputProcessor(fesys_controller, info, TECPLOT_OUTPUT_PROCESSOR::num()),
mesh(NULL),
mesh_data(NULL),
dummy_elem_set(NULL)
{
	
}




TecPlotOutputProcessor::~TecPlotOutputProcessor()
{
	
}



void TecPlotOutputProcessor::writeData()
{
//   // iterate over the disciplines and write the mesh to 
//   // each discipline's file
  
//   std::vector<Discipline::DisciplineInfo*>& disciplines = 
//   this->fesystem_controller.analysis_case->getAnalysisDisciplineInfos();
  
//   unsigned int mesh_ID =0;
  
//   std::vector<std::string> dof_names;

//   std::vector<Discipline::DisciplineInfo*>::const_iterator it, end;
//   it = disciplines.begin();
//   end = disciplines.end();
  
//   for (; it !=end ; it++)
//     {
//       mesh_ID = (*it)->getMeshID();
    
//       // get the mesh and meshdata from mesh list
//       this->mesh = 
// 	this->fesystem_controller.mesh_list->getMeshFromID(mesh_ID);
//       this->mesh_data = 
// 	this->fesystem_controller.mesh_list->getMeshDataFromID(mesh_ID);
//       this->dummy_elem_set = 
// 	this->fesystem_controller.mesh_list->getDummyElemSetFromID(mesh_ID);
      
//       // get the name of the file to be written to
//       std::string file_name = 
// 	Discipline::AnalysisDisciplineEnum::enumName((*it)->getDisciplineEnumID());
      
//       std::string name = this->output_info.getOutputFileName();
      
//       file_name += name;
      
//       // open the file, and call the functions to write the data to this stream
//       std::fstream output_file;
//       output_file.open(file_name.c_str(), std::fstream::out);
      
//       switch((*it)->getDisciplineEnumID())
// 	{
// 	case THERMAL_DISCIPLINE_ENUM_ID:
// 	  {
// 	    dof_names.clear();
// 	    dof_names.push_back("Temperature");
// 	  }
// 	  break;
	  
// 	case STRUCTURAL_DISCIPLINE_ENUM_ID:
// 	  {
// 	    dof_names.clear();
// 	    dof_names.push_back("u");
// 	    dof_names.push_back("v");
// 	    dof_names.push_back("w");
// 	    dof_names.push_back("theta_x");
// 	    dof_names.push_back("theta_y");
// 	    dof_names.push_back("theta_z");
// 	  }
// 	  break;
	  
// 	default:
// 	  Assert(false, ExcInternalError());
// 	}

//       // now write the nodal solution
//       this->writeNodalData(output_file, (*it)->getDisciplineEnumID(),
// 			   dof_names);
      
//       output_file.close();
//     }
}



void TecPlotOutputProcessor::prepareElemTypeMap()
{
//   // first of all, clear the map
//   this->elem_vec_map.clear();

//   // next, iterate over all the elems, and based on the elem type, 
//   // insert them in the vector. Finally, insert the vectors in the 
//   // map
//   // get the iterator for the elems
//   MeshBase::const_element_iterator elem_it = this->mesh->elements_begin();
//   MeshBase::const_element_iterator elem_end = this->mesh->elements_end();
  
//   std::set<Elem*>::const_iterator dummy_elem_end =
//   this->dummy_elem_set->end();
  
//   std::vector<Elem*> *elem_vec=NULL;
	
//   for (; elem_it != elem_end; elem_it++)
//     {
//       if (this->dummy_elem_set->find(*elem_it) != dummy_elem_end)
// 	continue;
      
      
//       // insert this elem type in the map if it does not already exist
//       if (this->elem_vec_map.count((**elem_it).type()) == 0)
// 	{
// 	  elem_vec = 
// 	    &(this->elem_vec_map.insert(std::map<ElemType, std::vector<Elem*> >::value_type
// 					((**elem_it).type(), std::vector<Elem*>())).first->second);
// 	}
//       else
// 	elem_vec = &(this->elem_vec_map.find((**elem_it).type())->second);
      
//       elem_vec->push_back(*elem_it);
//     }
  
	
//   for (; elem_it != elem_end; elem_it++)
//     {
//     Elem* elem = *elem_it;
    
//     if (this->dummy_elem_set->find(elem) != dummy_elem_end)
//       continue;
    
//     switch (elem->type())
//       {
//       case EDGE2:
//         {
//           edge2_elems.push_back(elem);
//         }
//         break;
        
//       case TRI3:
//         {
//           tri3_elems.push_back(elem);
//         }	
//         break;
        
//       case QUAD4:
//         {
//           quad4_elems.push_back(elem);
//         }	
//         break;
        
//       case HEX8:
//       case PRISM6:
//         {
//           hex8_elems.push_back(elem);
//         }	
//         break;
				
//       default:
//         Assert(false, ExcInternalError());
//         break;
//       }
//     }
  
//   // next, based on whether the elem vec is non-zero or not, insert them to the
//   // map
//   if (edge2_elems.size() > 0)
//     this->elem_vec_map.insert(std::map<ElemType, std::vector<Elem*> >::value_type
//                               (EDGE2, edge2_elems));
//   if (quad4_elems.size() > 0)
//     this->elem_vec_map.insert(std::map<ElemType, std::vector<Elem*> >::value_type
//                               (QUAD4, quad4_elems));
//   if (tri3_elems.size() > 0)
//     this->elem_vec_map.insert(std::map<ElemType, std::vector<Elem*> >::value_type
//                               (TRI3, tri3_elems));
//   if (hex8_elems.size() > 0)
//     this->elem_vec_map.insert(std::map<ElemType, std::vector<Elem*> >::value_type
//                               (HEX8, hex8_elems));
  
}






void TecPlotOutputProcessor::writeNodalData(std::ostream& output,
					    const unsigned int discipline,
					    const std::vector<std::string>& dof_names)
{
//   this->prepareElemTypeMap();

//   // get the solution data infos from all the solutions
  
  


//   // this is the vector in which the data will be obtained and then written from
//   std::auto_ptr<NumericVector<double> > vector
//   (NumericVector<double>::build().release()); 
  
  
//   // control parameters
//   output  << "TITLE = \" FESystem Output\" " << "\n"
//     << "VARIABLES = \"X\",  \"Y\",  \"Z\"";
  
//   // iterate over each solution, and prepare a list of nodal and element centered 
//   // variables for the same.
//   const std::vector<Solution::SolutionBase*> sols = 
//     this->fesystem_controller.analysis_case->getSolutions();

//   std::vector<Solution::SolutionBase*>::const_iterator sol_it, sol_end;
//   sol_it = sols.begin();
//   sol_end = sols.end();
  
//   std::vector<Solution::SolNameInfo>::const_iterator sol_name_it, sol_name_end;
//   std::vector<unsigned int>::const_iterator lc_it, lc_end;

//   unsigned int n_nodal_variables = 0;
  
//   // iterate over all the solutions and write the nodal solution names for the 
//   // load cases and DV sensitivities.
//   for ( ; sol_it != sol_end; sol_it++)
//     {
//     // if the solution does not have this discipline, nothing to be written.
//     if (! (*sol_it)->checkParticipatingDiscipline(Discipline::THERMAL_DISCIPLINE::num()))
//       continue;

//     // else, get the solution name and load cases from this solution
//     const std::vector<Solution::SolNameInfo >& sol_names = 
//       (*sol_it)->getNodalSolutionNames(Discipline::THERMAL_DISCIPLINE::num());
    
//     // iterate over all the solution names and write them to the output
//     sol_name_it = sol_names.begin();
//     sol_name_end = sol_names.end();
    
//     for (; sol_name_it != sol_name_end; sol_name_it++)
//       {
//       // get the load cases associated with this solution and write
//       // the names
//       lc_it = sol_name_it->load_cases.begin();
//       lc_end = sol_name_it->load_cases.end();
      
//       for ( ; lc_it != lc_end; lc_it++)
//         {
//         output << ", \"" << sol_name_it->sol_name << "_lc_" << *lc_it <<  "\"";
//         n_nodal_variables++;
//         }
      
//       }
//     }

//   output << "\n";

//   // having written all the variable names, write the zones and actual variable values
  
//   // also, get the iterators to the elem_type map
//   std::map<ElemType, std::vector<Elem*> >::const_iterator elem_map_it, elem_map_end;
//   elem_map_it = this->elem_vec_map.begin();
//   elem_map_end = this->elem_vec_map.end();
	
	
  
//   std::string zone_name, elem_name;
//   for (; elem_map_it != elem_map_end; elem_map_it++)
//     {
//     zone_name.clear();
//     elem_name.clear();
//     switch (elem_map_it->first)
//       {
//       case QUAD4:
//         {
//           zone_name = "QUAD ELEMS";
//           elem_name = "FEQUADRILATERAL";
//         }
//         break;
        
//       case EDGE2:
//         {
//           zone_name = "EDGE ELEMS";
//           elem_name = "FELINESEG";
//         }
//         break;
        
//       case TRI3:
//         {
//           zone_name = "TRI ELEMS";
//           elem_name = "FETRIANGLE";
//         }
//         break;
        
//       case HEX8:
//       case PRISM6:
//         {
//           zone_name = "HEX8 ELEMS";
//           elem_name = "FEBRICK";
//         }
//         break;
        
//       default:
//         abort();
//       }
    
    
//     if ( elem_map_it == this->elem_vec_map.begin())
//       {
//       output << "\n ZONE  T=\"" <<zone_name<< "\", N=" << this->mesh->n_nodes() << ",  E=" << elem_map_it->second.size() 
//       << ", DATAPACKING=BLOCK"<< ",  ZONETYPE="<< elem_name << std::endl;
//       // and write the node data
//       // node data
//       output.setf(std::ios::showpoint);
      
//       std::ostringstream node_x, node_y, node_z;
      
//       // iterate over the nodes and write the data
//       for (unsigned int node_it=0; node_it < this->mesh->n_nodes(); node_it++)
//         {
// 	      const Node* node = this->mesh_data->getNodeFromInternalID(node_it);
	      
// 	      // write the x,y,z data and then write the temperature data for each node
// 	      node_x << std::setprecision(8) << (*node)(0) << "  ";
// 	      node_y << std::setprecision(8) << (*node)(1) << "  ";
// 	      node_z << std::setprecision(8) << (*node)(2) << "  ";
//         }
      
//       node_x << std::endl; 
//       node_y << std::endl;
//       node_z << std::endl;
      
//       output << node_x.str();
//       output << node_y.str();
//       output << node_z.str();

//       // get the iterators to the solutions
//       sol_it = sols.begin();
//       sol_end = sols.end();
      
//       // iterate over all the solutions and write the nodal solutions for the 
//       // load cases
//       for ( ; sol_it != sol_end; sol_it++)
//         {
//         // if the solution does not have this discipline, nothing to be written.
//         if (! (*sol_it)->checkParticipatingDiscipline(Discipline::THERMAL_DISCIPLINE::num()))
//           continue;
        
//         // else, get the solution name and load cases from this solution
//         const std::vector<Solution::SolNameInfo >& sol_names = 
//           (*sol_it)->getNodalSolutionNames(Discipline::THERMAL_DISCIPLINE::num());
        
//         // iterate over all the solution names and write them to the output
//         sol_name_it = sol_names.begin();
//         sol_name_end = sol_names.end();
        
//         for (; sol_name_it != sol_name_end; sol_name_it++)
//           {
//           // get the load cases associated with this solution and write
//           // the names
//           lc_it = sol_name_it->load_cases.begin();
//           lc_end = sol_name_it->load_cases.end();
          
//           unsigned int dof=0;
//           for ( ; lc_it != lc_end; lc_it++)
//             {
//             // get the solution name in the database, fill the vector and write the variable
//             this->fesystem_controller.global_data_storage->fillVector(*lc_it,
//                                                                       sol_name_it->database_name,
//                                                                       *(vector.get()));
            
//             for (unsigned int node_it=0; node_it < this->mesh->n_nodes(); node_it++)
//               {
//               dof = this->mesh_data->getNodeFromInternalID(node_it)->dof_number(0,0,0);
//               output << std::setprecision(8) << (*vector.get())(dof) << "  ";
//               }
//             }
//           }
//         }
//       output << std::endl << std::endl;
//       }
//     else
//       {
//       output << "\n ZONE  T=\""<<zone_name<<"\", N=" << this->mesh->n_nodes() << ",  E=" << elem_map_it->second.size() 
//       << ", VARSHARELIST=([1-" 
//       << 3 + n_nodal_variables
//       <<"]=1)"<< ",  ZONETYPE="<< elem_name << std::endl;
//       }
	  
//     // element data
//     std::vector<Elem*>::const_iterator elem_it = elem_map_it->second.begin();
//     std::vector<Elem*>::const_iterator elem_end = elem_map_it->second.end();
    
//     unsigned int n_elem_nodes = 0;
//     for (; elem_it != elem_end; elem_it++)
//       {
//       n_elem_nodes = (*elem_it)->n_nodes();
//       const ElemType elemtype = (*elem_it)->type();
//       // for each node of the elem, find the node location in 
//       // the map, and write the connectivity data
//       switch (elemtype)
//         {
//         case PRISM6:
//           {
//             for (unsigned int i=0; i<n_elem_nodes; i++)
//               {
//               const Node* node = (*elem_it)->get_node(i);
//               output << (this->mesh_data->getInternalIDFromNode(node)+1) << "  " ;
//               if (i == 2 || i == 5)
//                 output << (this->mesh_data->getInternalIDFromNode(node)+1) << "  " ;
//               }
//           }
//           break;
          
//         default:
//           {
//             for (unsigned int i=0; i<n_elem_nodes; i++)
//               {
//               const Node* node = (*elem_it)->get_node(i);
//               output << (this->mesh_data->getInternalIDFromNode(node)+1) << "  " ;
//               }
//           }
//           break;
//         }
//       }
//     output << std::endl << std::endl << std::endl;
//     }
}



void TecPlotOutputProcessor::writeStructuralData(std::ostream& output)
{
//   this->prepareElemTypeMap();
  
// //  PostProcessQtyDatabase& database = 
// //    (*this->fesystem_controller.post_process_qty_database.get());
  
//   // this is the vector in which the data will be obtained and then written from
//   std::auto_ptr<NumericVector<double> > vector
//     (NumericVector<double>::build().release()); 
  
//   // also, get the iterators to the elem_type map
//   std::map<ElemType, std::vector<Elem*> >::const_iterator elem_map_it, elem_map_end;
//   elem_map_it = this->elem_vec_map.begin();
//   elem_map_end = this->elem_vec_map.end();
	
//   // control parameters
//   output  << "TITLE = \" FESystem Output\" " << "\n"
//     << "VARIABLES = \"X\",  \"Y\",  \"Z\"";

//   // iterate over each solution, and prepare a list of nodal and element centered 
//   // variables for the same.
//   const std::vector<Solution::SolutionBase*> sols = 
//     this->fesystem_controller.analysis_case->getSolutions();
  
//   std::vector<Solution::SolutionBase*>::const_iterator sol_it, sol_end;
//   sol_it = sols.begin();
//   sol_end = sols.end();
  
//   std::vector<Solution::SolNameInfo>::const_iterator sol_name_it, sol_name_end;
//   std::vector<unsigned int>::const_iterator lc_it, lc_end;
  
//   unsigned int n_nodal_variables = 0;
  
//   // iterate over all the solutions and write the nodal solution names for the 
//   // load cases and DV sensitivities.
//   for ( ; sol_it != sol_end; sol_it++)
//     {
//     // if the solution does not have this discipline, nothing to be written.
//     if (! (*sol_it)->checkParticipatingDiscipline(Discipline::STRUCTURAL_DISCIPLINE::num()))
//       continue;
    
//     // else, get the solution name and load cases from this solution
//     const std::vector<Solution::SolNameInfo >& sol_names = 
//       (*sol_it)->getNodalSolutionNames(Discipline::STRUCTURAL_DISCIPLINE::num());
    
//     // iterate over all the solution names and write them to the output
//     sol_name_it = sol_names.begin();
//     sol_name_end = sol_names.end();
    
//     for (; sol_name_it != sol_name_end; sol_name_it++)
//       {
//       // get the load cases associated with this solution and write
//       // the names
//       lc_it = sol_name_it->load_cases.begin();
//       lc_end = sol_name_it->load_cases.end();
      
//       for ( ; lc_it != lc_end; lc_it++)
//         {
//         output << ", \"" << sol_name_it->sol_name << "_u_lc_" << *lc_it <<  "\"";
//         output << ", \"" << sol_name_it->sol_name << "_v_lc_" << *lc_it <<  "\"";
//         output << ", \"" << sol_name_it->sol_name << "_w_lc_" << *lc_it <<  "\"";
//         output << ", \"" << sol_name_it->sol_name << "_theta_x_lc_" << *lc_it <<  "\"";
//         output << ", \"" << sol_name_it->sol_name << "_theta_y_lc_" << *lc_it <<  "\"";
//         output << ", \"" << sol_name_it->sol_name << "_theta_z_lc_" << *lc_it <<  "\"";
//         n_nodal_variables += 6;
//         }
      
//       }
//     }
  
// //  for (unsigned int i=0; i<load_cases.size(); i++)
// //    {
// //    output << ", \"epsilon_xx_lc_" << load_cases[i] <<  "\"";
// //    output << ", \"epsilon_yy_lc_" << load_cases[i] <<  "\"";
// //    output << ", \"epsilon_xy_lc_" << load_cases[i] <<  "\"";
// //    output << ", \"sigma_xx_lc_" << load_cases[i] <<  "\"";
// //    output << ", \"sigma_yy_lc_" << load_cases[i] <<  "\"";
// //    output << ", \"sigma_xy_lc_" << load_cases[i] <<  "\"";
// //    output << ", \"strain_energy_density_lc_" << load_cases[i] <<  "\"";
// //    output << ", \"VonMises_lc_" << load_cases[i] <<  "\"";
// //    }
//   output << "\n";
	
  
//   std::string zone_name, elem_name;
//   for (; elem_map_it != elem_map_end; elem_map_it++)
//     {
//     zone_name.clear();
//     elem_name.clear();
//     switch (elem_map_it->first)
//       {
//       case QUAD4:
//         {
//           zone_name = "QUAD ELEMS";
//           elem_name = "FEQUADRILATERAL";
//         }
//         break;
        
//       case EDGE2:
//         {
//           zone_name = "EDGE ELEMS";
//           elem_name = "FELINESEG";
//         }
//         break;
        
//       case TRI3:
//         {
//           zone_name = "TRI ELEMS";
//           elem_name = "FETRIANGLE";
//         }
//         break;
        
//       case HEX8:
//       case PRISM6:
//         {
//           zone_name = "HEX8 ELEMS";
//           elem_name = "FEBRICK";
//         }
//         break;
        
//       default:
//         abort();
//       }
    
//     if ( elem_map_it == this->elem_vec_map.begin())
//       {
// //      output << "\n ZONE  T=\""<<zone_name<<"\", N=" << this->mesh->n_nodes() << ",  E=" << elem_map_it->second.size() 
// //      << ", DATAPACKING=BLOCK"<< ", VARLOCATION=([1-"
// //      << 3+6*n_load_cases<<"]=NODAL,["
// //      << 3+6*n_load_cases+1<<"-"
// //      << 3+(6+8)*n_load_cases<<"]=CELLCENTERED)" 
// //      << ", ZONETYPE="<< elem_name << std::endl;

//       output << "\n ZONE  T=\""<<zone_name<<"\", N=" << this->mesh->n_nodes() << ",  E=" << elem_map_it->second.size() 
//       << ", DATAPACKING=BLOCK"<< ", VARLOCATION=([1-"
//       << 3+n_nodal_variables<<"]=NODAL)" 
//       << ", ZONETYPE="<< elem_name << std::endl;
      
//       // and write the node data
//       // node data
//       output.setf(std::ios::showpoint);
      
//       std::ostringstream node_x, node_y, node_z;
      
//       // iterate over the nodes and write the data
//       for (unsigned int node_it=0; node_it < this->mesh->n_nodes(); node_it++)
//         {
// 	      const Node* node = this->mesh_data->getNodeFromInternalID(node_it);
	      
// 	      // write the x,y,z data and then write the temperature data for each node
// 	      node_x << std::setprecision(8) << (*node)(0) << "  ";
// 	      node_y << std::setprecision(8) << (*node)(1) << "  ";
// 	      node_z << std::setprecision(8) << (*node)(2) << "  ";
//         }
      
//       node_x << std::endl; 
//       node_y << std::endl;
//       node_z << std::endl;
      
//       output << node_x.str();
//       output << node_y.str();
//       output << node_z.str();
      
//       // get the iterators to the solutions
//       sol_it = sols.begin();
//       sol_end = sols.end();
      
//       // iterate over all the solutions and write the nodal solutions for the 
//       // load cases
//       for ( ; sol_it != sol_end; sol_it++)
//         {
//         // if the solution does not have this discipline, nothing to be written.
//         if (! (*sol_it)->checkParticipatingDiscipline(Discipline::STRUCTURAL_DISCIPLINE::num()))
//           continue;
        
//         // else, get the solution name and load cases from this solution
//         const std::vector<Solution::SolNameInfo >& sol_names = 
//           (*sol_it)->getNodalSolutionNames(Discipline::STRUCTURAL_DISCIPLINE::num());
        
//         // iterate over all the solution names and write them to the output
//         sol_name_it = sol_names.begin();
//         sol_name_end = sol_names.end();
        
//         for (; sol_name_it != sol_name_end; sol_name_it++)
//           {
//           // get the load cases associated with this solution and write
//           // the names
//           lc_it = sol_name_it->load_cases.begin();
//           lc_end = sol_name_it->load_cases.end();
          
//           unsigned int dof=0;
//           for ( ; lc_it != lc_end; lc_it++)
//             {
//             std::ostringstream u_var, v_var, w_var, theta_x_var, theta_y_var, theta_z_var;

//             // get the solution name in the database, fill the vector and write the variable
//             this->fesystem_controller.global_data_storage->fillVector(*lc_it,
//                                                                       sol_name_it->database_name,
//                                                                       *(vector.get()));
            
//             for (unsigned int node_it=0; node_it < this->mesh->n_nodes(); node_it++)
//               {
//               for (unsigned int i=0; i<6; i++)
//                 {
//                 dof = this->mesh_data->getNodeFromInternalID(node_it)->dof_number(0,i,0);
//                 switch (i)
//                   {
//                   case 0:
//                     u_var << std::setprecision(8) << (*vector.get())(dof) << "  ";
//                     break;
//                   case 1:
//                     v_var << std::setprecision(8) << (*vector.get())(dof) << "  ";
//                     break;
//                   case 2:
//                     w_var << std::setprecision(8) << (*vector.get())(dof) << "  ";
//                     break;
//                   case 3:
//                     theta_x_var << std::setprecision(8) << (*vector.get())(dof) << "  ";
//                     break;
//                   case 4:
//                     theta_y_var << std::setprecision(8) << (*vector.get())(dof) << "  ";
//                     break;
//                   case 5:
//                     theta_z_var << std::setprecision(8) << (*vector.get())(dof) << "  ";
//                     break;
//                   default:
//                     abort();
//                     break;
//                   }
//                 }
//               }
            
//             output << u_var.str() << std::endl;
//             output << v_var.str() << std::endl;
//             output << w_var.str() << std::endl;
//             output << theta_x_var.str() << std::endl;
//             output << theta_y_var.str() << std::endl;
//             output << theta_z_var.str() << std::endl;
//             }
//           }
//         }
//       output << std::endl << std::endl;
//       }
//     else
//       {
// //      output << "\n ZONE  T=\""<< zone_name << "\", N=" << this->mesh->n_nodes() << ",  E=" << elem_map_it->second.size() 
// //      << ", DATAPACKING=BLOCK"<< ", VARLOCATION=([1-"
// //      << 3+6*n_load_cases<<"]=NODAL,["
// //      << 3+6*n_load_cases+1<<"-"
// //      << 3+(6+8)*n_load_cases<<"]=CELLCENTERED)" 
// //      << ", DATAPACKING=BLOCK" 
// //      <<  ", VARSHARELIST=([1-"
// //      << 3+6*n_load_cases << "]=1)"<< ",  ZONETYPE=" << elem_name <<" \n";
//       output << "\n ZONE  T=\""<< zone_name << "\", N=" << this->mesh->n_nodes() << ",  E=" << elem_map_it->second.size() 
//       << ", DATAPACKING=BLOCK"<< ", VARLOCATION=([1-"
//       << 3+n_nodal_variables<<"]=NODAL)" 
//       << ", DATAPACKING=BLOCK" 
//       <<  ", VARSHARELIST=([1-"
//       << 3+n_nodal_variables << "]=1)"<< ",  ZONETYPE=" << elem_name <<" \n";
//       }
	  
// //    // write the element dependent variables
// //    unsigned int elem_ID;
// //    lc_it = load_cases.begin();
// //    for ( ; lc_it != lc_end; lc_it++)
// //      {
// //      std::ostringstream lc_num;
// //      lc_num << *lc_it;
// //	    
// //      std::ostringstream 
// //        stream_vector_0, stream_vector_1, stream_vector_2, stream_vector_3, 
// //        stream_vector_4, stream_vector_5, stream_vector_6, stream_vector_7;
// //      
// //      // get the iterator for the elems
// //      std::vector<Elem*>::const_iterator elem_it = elem_map_it->second.begin();
// //      std::vector<Elem*>::const_iterator elem_end = elem_map_it->second.end();
// //      
// //      std::set<Elem*>::const_iterator dummy_elem_end =
// //        this->dummy_elem_set->end();
// //	    
// //      for (; elem_it != elem_end; elem_it++)
// //        {
// //	      Elem* elem = *elem_it;
// //        
// //	      if (this->dummy_elem_set->find(elem) != dummy_elem_end)
// //          continue;
// //        
// //        
// //	      elem_ID = this->mesh_data->getForeignIDFromElem(elem);
// //        
// //	      ElemPostProcessQty& qty =  
// //          dynamic_cast<ElemPostProcessQty&>(database.getElementPostProcessQty
// //                                            (Discipline::STRUCTURAL_DISCIPLINE::num(),
// //                                             elem_ID));
// //        
// //	      // get the strain and stress tensor
// //	      TensorValue<double>& strain = qty.getStrainTensor(*lc_it);
// //	      TensorValue<double>& stress = qty.getStressTensor(*lc_it);
// //        
// //	      // write the values to the streams
// //	      stream_vector_0 << strain(0,0) << "  ";
// //	      stream_vector_1 << strain(1,1) << "  ";
// //	      stream_vector_2 << strain(0,1) << "  ";
// //        
// //	      stream_vector_3 << stress(0,0) << "  ";
// //	      stream_vector_4 << stress(1,1) << "  ";
// //	      stream_vector_5 << stress(0,1) << "  ";
// //        
// //	      stream_vector_6 << qty.getStrainEnergyDensity(*lc_it) << "  ";
// //	      stream_vector_7 << qty.getVonMisesStress(*lc_it) << "  " ;
// //        }	
// //	    
// //      // now write the variables to the output
// //      output << stream_vector_0.str() << std::endl;
// //      output << stream_vector_1.str() << std::endl;
// //      output << stream_vector_2.str() << std::endl;
// //      output << stream_vector_3.str() << std::endl;
// //      output << stream_vector_4.str() << std::endl;
// //      output << stream_vector_5.str() << std::endl;
// //      output << stream_vector_6.str() << std::endl;
// //      output << stream_vector_7.str() << std::endl;
// //      }
    
// 	  // element data
// 	  std::vector<Elem*>::const_iterator elem_it = elem_map_it->second.begin();
// 	  std::vector<Elem*>::const_iterator elem_end = elem_map_it->second.end();
	  
// 	  unsigned int n_elem_nodes = 0;
// 	  for (; elem_it != elem_end; elem_it++)
// 	    {
//       n_elem_nodes = (*elem_it)->n_nodes();
//       const ElemType elemtype = (*elem_it)->type();
//       // for each node of the elem, find the node location in 
//       // the map, and write the connectivity data
//       switch (elemtype)
//         {
//         case PRISM6:
//           {
//             for (unsigned int i=0; i<n_elem_nodes; i++)
//               {
//               const Node* node = (*elem_it)->get_node(i);
//               output << (this->mesh_data->getInternalIDFromNode(node) + 1) << "  " ;
//               if (i == 2 || i == 5)
//                 output << (this->mesh_data->getInternalIDFromNode(node) + 1) << "  " ;
//               }
//           }
//           break;
          
//         default:
//           {
//             for (unsigned int i=0; i<n_elem_nodes; i++)
//               {
//               const Node* node = (*elem_it)->get_node(i);
//               output << (this->mesh_data->getInternalIDFromNode(node) + 1) << "  " ;
//               }
//           }
//           break;
//         }
// 	    }
// 	  output << std::endl << std::endl << std::endl;
//     }
}




